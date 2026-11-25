#pragma once

#include <volt/core/particle_property.h>
#include <volt/core/progress.h>
#include <volt/core/simulation_cell.h>
#include <volt/math/matrix3.h>
#include <volt/math/quaternion.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/analysis/ptm.h>
#include <volt/analysis/nearest_neighbor_finder.h>

#include <ptm_functions.h>
#include <boost/sort/sort.hpp>
#include <boost/intrusive/rbtree_algorithms.hpp>
#include <boost/heap/priority_queue.hpp>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/partitioner.h>
#include <tbb/parallel_sort.h>

#include <vector>
#include <unordered_set>
#include <queue>
#include <cmath>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <algorithm>

namespace Volt{

class GrainSegmentationLog;

enum class MergeAlgorithm{
    GraphClusteringAutomatic,
    GraphClusteringManual,
    MinimumSpanningTree
};

inline const char* mergeAlgorithmName(MergeAlgorithm algorithm){
    switch(algorithm){
        case MergeAlgorithm::GraphClusteringAutomatic:
            return "GraphClusteringAutomatic";
        case MergeAlgorithm::GraphClusteringManual:
            return "GraphClusteringManual";
        case MergeAlgorithm::MinimumSpanningTree:
            return "MinimumSpanningTree";
    }
    return "Unknown";
}

union NodeUnion{
    size_t opposite;
    size_t size;
};

struct HalfEdge{
    HalfEdge* _parent;
    HalfEdge* _left;
    HalfEdge* _right;
    int _color;

    NodeUnion data;
    double weight;
};

inline double calculateMedian(std::vector<double>& data){
    size_t n = data.size();
    std::sort(data.begin(), data.end());
    double median = data[n / 2];
    if(n % 2 == 0){
        median += data[n / 2 - 1];
        median /= 2;
    }
    return median;
}

inline void weightedLinearRegression(
    std::vector<double>& weights,
    std::vector<double>& xs,
    std::vector<double>& ys,
    double& gradient,
    double& intercept
){
    double wsum = 0;
    for(auto w : weights){
        wsum += w;
    }

    for(size_t i = 0; i < weights.size(); i++){
        weights[i] /= wsum;
    }

    double xmean = 0;
    double ymean = 0;
    for(size_t i = 0; i < weights.size(); i++){
        xmean += weights[i] * xs[i];
        ymean += weights[i] * ys[i];
    }

    double sum_xx = 0;
    double sum_xy = 0;
    for(size_t i = 0; i < weights.size(); i++){
        sum_xx += weights[i] * (xs[i] - xmean) * (xs[i] - xmean);
        sum_xy += weights[i] * (xs[i] - xmean) * (ys[i] - ymean);
    }

    gradient = sum_xy / sum_xx;
    intercept = ymean - gradient * xmean;
}

inline std::vector<double> leastAbsoluteDeviations(
    std::vector<double>& weights, 
    std::vector<double>& xs, 
    std::vector<double>& ys, 
    double& gradient, 
    double& intercept
){
    std::vector<double> residuals(weights.size());
    std::vector<double> w(weights);

    for(int it = 0; it < 100; it++){
        weightedLinearRegression(w, xs, ys, gradient, intercept);
        
        for(size_t i = 0; i < xs.size(); i++){
            double prediction = gradient * xs[i] + intercept;
            double r = std::abs(ys[i] - prediction);
            residuals[i] = r;
            w[i] = weights[i] / std::max(1E-4, r);
        }
    }

    return residuals;
}  

struct rbtree_node_traits{
    typedef HalfEdge node;
    typedef HalfEdge* node_ptr;
    typedef const HalfEdge* const_node_ptr;
    typedef int color;

    static node_ptr get_parent(const_node_ptr n){
        return n->_parent;
    }

    static void set_parent(node_ptr n, node_ptr p){
        n->_parent = p;
    }

    static node_ptr get_left(const_node_ptr n){
        return n->_left;
    }
    
    static void set_left(node_ptr n, node_ptr l){
        n->_left = l;
    }

    static node_ptr get_right(const_node_ptr n){
        return n->_right;
    }

    static void set_right(node_ptr n, node_ptr r){
        n->_right = r;
    }

    static color get_color(const_node_ptr n){
        return n->_color;
    }

    static void set_color(node_ptr n, color c){
        n->_color = c;
    }

    static color black(){
        return color(0);
    }

    static color red(){
        return color(1);
    }
};

struct node_ptr_compare{
    bool operator()(const HalfEdge* a, const HalfEdge* b) const{
        return a->data.opposite < b->data.opposite;
    }
};

typedef boost::intrusive::rbtree_algorithms<rbtree_node_traits> algo;

class DisjointSet{
public:
    DisjointSet(size_t n){
        parents.resize(n);
        sizes.resize(n);
        clear();
    }

    void clear(){
        std::iota(parents.begin(), parents.end(), (size_t) 0);
        std::fill(sizes.begin(), sizes.end(), 1);
    }

    size_t find(size_t index){
        size_t x = parents[index];
        while(x != parents[x]){
            parents[x] = parents[parents[x]];
            x = parents[x];
        }
        parents[index] = x;
        return x;
    }

    size_t merge(size_t index1, size_t index2){
        size_t parentA = find(index1);
        size_t parentB = find(index2);
        if(parentA == parentB) return parentA;

        if(sizes[parentA] < sizes[parentB]){
            parents[parentA] = parentB;
            sizes[parentB] += sizes[parentA];
            return parentB;
        }else{
            parents[parentB] = parentA;
            sizes[parentA] += sizes[parentB];
            return parentA;
        }
    }

    size_t nodesize(size_t index) const{
        return sizes[index];
    }
    
private:
    std::vector<size_t> parents;
    std::vector<size_t> sizes;
};

class GrainSegmentationEngine1{
public:
    class Graph{
    public:
        std::vector<double> wnode;
        std::vector<HalfEdge> header;
        std::vector<HalfEdge> edgeBuffer;
        size_t edgeCount = 0;
        std::unordered_set<size_t> activeNodes;

        Graph(size_t numNodes, size_t numEdges){
            wnode.assign(numNodes, 0.0);
            header.resize(numNodes);
            for(size_t i = 0; i < numNodes; ++i){
                algo::init_header(&header[i]);
                header[i].data.size = 0;
            }
            edgeBuffer.resize(2 * numEdges);
        }

        size_t num_nodes() const{
            return activeNodes.size();
        }

        size_t next_node() const{
            return *activeNodes.begin();
        } 

        std::tuple<double, size_t> nearestNeighbor(size_t a) const{
            double dmin = std::numeric_limits<double>::infinity();
            size_t vmin = std::numeric_limits<size_t>::max();

            HalfEdge* e = header[a]._left;
            for(size_t it = 0; it < header[a].data.size; ++it){
                size_t v = e->data.opposite;
                double w = e->weight;
                e = algo::next_node(e);

                if(v == a){
                    throw std::runtime_error("Graph has self loops");
                }

                double d = wnode[v] / std::max(w, 1e-300);
                if(d < dmin || (d == dmin && v < vmin)){
                    dmin = d; vmin = v;
                }
            }

            return std::make_tuple(dmin * wnode[a], vmin);
        }

        void add_edge(size_t u, size_t v, double w){
            size_t nodes[2] = {u,v};
            for(size_t idx : nodes){
                if(header[idx].data.size == 0){
                    activeNodes.insert(idx);
                }

                wnode[idx] += w;
            }

            HalfEdge* e = &edgeBuffer[edgeCount++];
            e->data.opposite = v; e->weight = w;
            insert_halfedge(&header[u], e);

            e = &edgeBuffer[edgeCount++];
            e->data.opposite = u; e->weight = w;
            insert_halfedge(&header[v], e);
        }

        void remove_node(size_t u){
            activeNodes.erase(u);
        }

        size_t contract_edge(size_t a, size_t b){
            if(header[b].data.size > header[a].data.size){
                std::swap(a,b);
            }

            algo::unlink(find(&header[b], a));
            algo::unlink(find(&header[a], b));
            header[a].data.size--;
            header[b].data.size--;

            HalfEdge* edge = header[b]._left;
            while(header[b].data.size != 0){
                size_t v = edge->data.opposite;
                double w = edge->weight;
                HalfEdge* next = algo::next_node(edge);
                algo::unlink(edge);
                header[b].data.size--;

                HalfEdge* opposite = find(&header[v], b);
                algo::unlink(opposite);
                header[v].data.size--;

                HalfEdge* temp = find(&header[a], v);
                if(temp){
                    temp->weight += w;
                    find(&header[v], a)->weight += w;
                }else{
                    edge->data.opposite = v; edge->weight = w;
                    insert_halfedge(&header[a], edge);

                    opposite->data.opposite = a; opposite->weight = w;
                    insert_halfedge(&header[v], opposite);
                }

                edge = next;
            }

            remove_node(b);
            wnode[a] += wnode[b];
            return a;
        }

    private:
        static HalfEdge* find(HalfEdge* header, size_t index){
            HalfEdge key;
            key.data.opposite = index;
            HalfEdge* result = algo::find(header, &key, node_ptr_compare());
            return result == header ? nullptr : result;
        }

        static void insert_halfedge(HalfEdge* header, HalfEdge* edge){
            algo::insert_equal_upper_bound(header, edge, node_ptr_compare());
            header->data.size++;
        }
    };

    static constexpr int MAX_DISORDERED_NEIGHBORS = 8;

    struct NeighborBond{
        size_t a;
        size_t b;
        double disorientation;
        double length;
    };

    struct DendrogramNode{
        DendrogramNode() = default;
        DendrogramNode(size_t _a, size_t _b, double _distance, double _disorientation, size_t _size, const Quaternion& _orientation)
            : a(_a)
            , b(_b)
            , distance(_distance)
            , disorientation(_disorientation)
            , size(_size)
            , merge_size(0.0)
            , orientation(_orientation){}
        size_t a = 0;
        size_t b = 0;
        double distance = 0.0;
        double disorientation = 0.0;
        size_t size = 0;
        double merge_size = 0.0;
        Quaternion orientation;
    };

    // The bonds between each atom and its nearest neighbors, before their disorientations are known.
    struct NeighborBonds{
        const std::vector<NeighborBond>& bonds;
        size_t atoms;
    };

    // The atoms of the minority phase of each coherent interface rotated onto the parent phase, and the parent of
    // the FCC/HCP pair.
    struct CoherentInterfaces{
        size_t rotatedAtoms;
        StructureType fccParent;
    };

    // The bonds once each has its disorientation, sorted by it.
    struct BondDisorientations{
        const std::vector<NeighborBond>& bonds;
    };

    class InterfaceHandler{
    public:
        explicit InterfaceHandler(const std::shared_ptr<ParticleProperty>& structures){
            int counts[(int) StructureType::NUM_STRUCTURE_TYPES] = {0};

            for(size_t i = 0; i < structures->size(); ++i){
                int t = structures->getInt(i);
                if(t >= 0 && t < (int) StructureType::NUM_STRUCTURE_TYPES){
                    counts[t]++;
                }
            }

            parent_fcc = counts[(int) StructureType::FCC] >= counts[(int) StructureType::HCP];
            parent_dcub = counts[(int) StructureType::CUBIC_DIAMOND] >= counts[(int) StructureType::HEX_DIAMOND];

            for(int k = 0;k < (int) StructureType::NUM_STRUCTURE_TYPES; ++k){
                target[k] = (StructureType) k;
            }

            if(parent_fcc){
                target[(int) StructureType::HCP] = StructureType::FCC;
            }else{
                target[(int) StructureType::FCC] = StructureType::HCP;
            }

            if(parent_dcub){
                target[(int) StructureType::HEX_DIAMOND] = StructureType::CUBIC_DIAMOND;
            }else{
                target[(int) StructureType::CUBIC_DIAMOND] = StructureType::HEX_DIAMOND;
            }
        }

        StructureType parent_phase(StructureType s) const{
            return target[(int) s];
        }

        bool reorder_bond(NeighborBond& bond, const std::vector<StructureType>& types) const{
            size_t a = bond.a;
            size_t b = bond.b;

            auto sa = types[a];
            auto sb = types[b];

            bool flipped = false;
            if(sa == StructureType::FCC && sb == StructureType::HCP){
                flipped |= !parent_fcc;
            }else if(sa == StructureType::HCP && sb == StructureType::FCC){
                flipped |=  parent_fcc;
            }else if(sa == StructureType::CUBIC_DIAMOND && sb == StructureType::HEX_DIAMOND){
                flipped |= !parent_dcub;
            }else if(sa == StructureType::HEX_DIAMOND && sb == StructureType::CUBIC_DIAMOND){
                flipped |=  parent_dcub;
            }else{
                return false;
            }

            if(flipped) std::swap(a,b);

            bond.a = a; 
            bond.b = b;

            return true;
        }

    private:
        bool parent_fcc  = true;
        bool parent_dcub = true;
        StructureType target[(int) StructureType::NUM_STRUCTURE_TYPES];
    };

    class Regressor{
    public:
        double gradient = 0;
        double intercept = 0;
        double mean_absolute_deviation = 0;

        std::vector<double> residuals;
        std::vector<double> xs;
        std::vector<double> ys;
        std::vector<double> weights;

        Regressor(std::vector<GrainSegmentationEngine1::DendrogramNode>& dendrogram){
            if(dendrogram.size() == 0){
                return;
            }

            for(auto node : dendrogram){
                weights.push_back(node.merge_size);
                xs.push_back(log(node.merge_size));
                ys.push_back(log(node.distance));
            }

            residuals = leastAbsoluteDeviations(weights, xs, ys, gradient, intercept);
            mean_absolute_deviation = calculateMedian(residuals); 
        }

        double calculate_threshold(std::vector<GrainSegmentationEngine1::DendrogramNode>& dendrogram, double cutoff){
            double threshold = 0;
            for(auto node : dendrogram) {
                double x = log(node.merge_size);
                double y = log(node.distance);

                double prediction = x * gradient + intercept;
                double residual = y - prediction;
                if (residual < cutoff * mean_absolute_deviation) {
                    threshold = std::max(threshold, y);
                }
            }

            return threshold;
        }
    };

    GrainSegmentationEngine1(
        std::shared_ptr<ParticleProperty> positions,
        std::shared_ptr<ParticleProperty> structures,
        std::shared_ptr<ParticleProperty> orientations,
        std::shared_ptr<ParticleProperty> correspondences,
        const SimulationCell* simCell,
        MergeAlgorithm algorithmType,
        bool handleCoherentInterfaces,
        bool outputBonds
    )
    : _algorithmType(algorithmType)
    , _handleBoundaries(handleCoherentInterfaces)
    , _numParticles(positions ? positions->size() : 0)
    , _positions(std::move(positions))
    , _structuresProperty(std::move(structures))
    , _orientationsProperty(std::move(orientations))
    , _correspondencesProperty(std::move(correspondences))
    , _simCell(*simCell)
    , _outputBonds(outputBonds)
    {
        _adjustedStructureTypes.resize(_numParticles, StructureType::OTHER);
        _adjustedOrientations.resize(_numParticles);

        for(size_t i=0;i<_numParticles;++i){
            _adjustedStructureTypes[i] = (StructureType)_structuresProperty->getInt(i);
            const double* q = _orientationsProperty->dataDouble() + 4 * i;
            _adjustedOrientations[i] = Quaternion(q[0], q[1], q[2], q[3]);
            _adjustedOrientations[i].normalize();
        }
    }

    void perform(GrainSegmentationLog& log);

    const std::vector<DendrogramNode>& dendrogram() const{
        return _dendrogram;
    }

    double suggestedMergingThreshold() const{
        return _suggestedMergingThreshold;
    }

    MergeAlgorithm algorithmType() const{
        return _algorithmType;
    }

    const std::shared_ptr<ParticleProperty>& structuresProperty() const{
        return _structuresProperty;
    }

    const std::shared_ptr<ParticleProperty>& orientationsProperty() const{
        return _orientationsProperty;
    }

    const std::vector<NeighborBond>& neighborBonds() const{
        return _neighborBonds;
    }

    // Whether the bond is an edge of the graph the merge sequence clusters.
    bool joinsMergeGraph(const NeighborBond& bond) const{
        return isCrystallineBond(bond) && bond.disorientation < _misorientationThresholdDeg;
    }

    static constexpr double misorientationThresholdDeg(){
        return _misorientationThresholdDeg;
    }

private:
    static inline double calculateGraphWeight(double theta_deg){
        if(theta_deg < 1e-5){
            theta_deg = 0.0;
        }

        return std::exp(-(1.0/3.0) * theta_deg * theta_deg);
    }

    static inline int desired_ptm_neighbor_count(StructureType st, int available){
        if(st == StructureType::OTHER){
            return std::min(available, MAX_DISORDERED_NEIGHBORS);
        }

        int ptmType = PTM::toPtmStructureType(st);
        int want = ptm_num_nbrs[ptmType];
        return std::min(available, want);
    } 

private:
    void createNeighborBonds(){
        ProgressSteps steps({
            {"Finding nearest neighbors", 1.0},
            {"Building neighbor bonds", 1.0},
        });
        PTM neighFinder;
        if(!neighFinder.prepare(_positions->constDataPoint3(), _numParticles, _simCell)){
            throw std::runtime_error("Error trying to prepare PTM neighbor finder.");
        }
        steps.next();

        using BaseQuery = NearestNeighborFinder::Query<PTM::MAX_INPUT_NEIGHBORS>;
        tbb::enumerable_thread_specific<BaseQuery> baseQueries([&]{
            return BaseQuery(neighFinder);
        });

        tbb::enumerable_thread_specific<std::vector<NeighborBond>> tlsBonds;

        ProgressReporter progress("  Building neighbor bonds", _numParticles);
        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, _numParticles, 1024), [&](const tbb::blocked_range<std::size_t>& r){
            auto& base = baseQueries.local();
            auto& out = tlsBonds.local();

            for(size_t i = r.begin(); i != r.end(); ++i){
                base.findNeighbors(i);
                const auto& res = base.results();
                int available = (int) res.size();

                StructureType st = (StructureType) _structuresProperty->getInt(i);
                int num = desired_ptm_neighbor_count(st, available);

                for(int j = 0; j < num; ++j){
                    size_t nb = res[j].index;
                    if(i < nb){
                        double length = std::sqrt(res[j].distanceSq);
                        out.push_back({ i, nb, std::numeric_limits<double>::infinity(), length });
                    }

                    const Vector3& neighborVector = res[j].delta;
                    for(size_t dim = 0; dim < 3; ++dim){
                        if(_simCell.hasPbc(dim)){
                            if(std::abs(_simCell.inverseMatrix().prodrow(neighborVector, dim)) >= 0.5 + EPSILON){
                                static const char* axes[3] = { "X", "Y", "Z" };
                                throw std::runtime_error(
                                    std::string("Simulation box is too short along cell vector ")
                                    + std::to_string(dim + 1) + " (" + axes[dim] + ") to perform "
                                    "grain segmentation. Please extend it first using the 'Replicate' modifier.");
                            }
                        }
                    }
                }
            }
            progress.advance(r.size());
        }, tbb::auto_partitioner{});

        tlsBonds.combine_each([&](const std::vector<NeighborBond>& v){
            _neighborBonds.insert(_neighborBonds.end(), v.begin(), v.end());
        });
    }

    bool interface_cubic_hex(NeighborBond& bond, const InterfaceHandler& iface, Quaternion& outRot){
        bond.disorientation = std::numeric_limits<double>::infinity();
        if(!iface.reorder_bond(bond, _adjustedStructureTypes)){
            return false;
        }

        size_t a = bond.a;
        size_t b = bond.b;

        const StructureType sa = _adjustedStructureTypes[a];
        const StructureType sb = _adjustedStructureTypes[b];

        Quaternion qrot;
        double mis = PTM::calculateInterfacialDisorientation(sa, sb, _adjustedOrientations[a], _adjustedOrientations[b], qrot);
        bond.disorientation = mis;
        outRot = qrot;
        return mis < _misorientationThresholdDeg;
    }

    CoherentInterfaces rotateInterfaceAtoms(){
        if(_neighborBonds.empty()) createNeighborBonds();

        InterfaceHandler iface(_structuresProperty);

        PTM neighFinder;
        if(!neighFinder.prepare(_positions->constDataPoint3(), _numParticles, _simCell)){
            throw std::runtime_error("Error trying to prepare PTM neighbor finder.");
        }

        using BaseQuery = NearestNeighborFinder::Query<PTM::MAX_INPUT_NEIGHBORS>;
        BaseQuery base(neighFinder);

        struct PQCmp{
            bool operator()(const NeighborBond& a, const NeighborBond& b) const{
                return a.disorientation > b.disorientation; 
            }
        };

        std::priority_queue<NeighborBond, std::vector<NeighborBond>, PQCmp> pq;
        size_t rotatedAtoms = 0;

        for(auto b : _neighborBonds){
            Quaternion rot;
            if(interface_cubic_hex(b, iface, rot)){
                pq.push(b);
            }
        }

        while(!pq.empty()){
            auto bond = pq.top(); 
            pq.pop();

            Quaternion rotated;
            if(!interface_cubic_hex(bond, iface, rotated)) continue;

            size_t idx = bond.b; 
            _adjustedStructureTypes[idx] = iface.parent_phase(_adjustedStructureTypes[idx]);
            _adjustedOrientations[idx]   = rotated;
            ++rotatedAtoms;

            base.findNeighbors(idx);
            const auto& res = base.results();
            int available = (int) res.size();
            int num = desired_ptm_neighbor_count(_adjustedStructureTypes[idx], available);

            for(int j = 0;j < num; ++j){
                size_t nb = res[j].index;
                NeighborBond b2{ idx, nb, 0.0, std::sqrt(res[j].distanceSq) };
                if(interface_cubic_hex(b2, iface, rotated)){
                    pq.push(b2);
                }
            }
        }
        return {rotatedAtoms, iface.parent_phase(StructureType::FCC)};
    }

    bool isCrystallineBond(const NeighborBond& b) const{
        auto a = _adjustedStructureTypes[b.a];
        auto c = _adjustedStructureTypes[b.b];

        if(a == StructureType::OTHER) return false;
        if(c == StructureType::OTHER) return false;
        if(a == c) return true;
        if(!_handleBoundaries) return false;

        if((a == StructureType::FCC && c == StructureType::HCP) || (a == StructureType::HCP && c == StructureType::FCC)) return true;
        if((a == StructureType::CUBIC_DIAMOND && c == StructureType::HEX_DIAMOND) || (a == StructureType::HEX_DIAMOND && c == StructureType::CUBIC_DIAMOND)) return true;
        return false;
    }

    void computeDisorientationAngles(){
        if(_neighborBonds.empty()) createNeighborBonds();
        const size_t N = _neighborBonds.size();

        InterfaceHandler iface(_structuresProperty);

        ProgressReporter progress("  Computing bond disorientations", N);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, N, 1024), [&](const tbb::blocked_range<size_t>& r){
            for(size_t i = r.begin(); i != r.end(); ++i){
                auto& b = _neighborBonds[i];
                auto sa = _adjustedStructureTypes[b.a];
                auto sb = _adjustedStructureTypes[b.b];

                if(sa == StructureType::OTHER || sb == StructureType::OTHER){
                    b.disorientation = std::numeric_limits<double>::infinity();
                    continue;
                }

                if(sa == sb){
                    b.disorientation = PTM::calculateDisorientation(sa, sb, _adjustedOrientations[b.a], _adjustedOrientations[b.b]);
                }else if(_handleBoundaries){
                    Quaternion dummy;
                    auto tmp = b;
                    if(interface_cubic_hex(tmp, iface, dummy)){
                        b.disorientation = tmp.disorientation;
                    }else{
                        b.disorientation = std::numeric_limits<double>::infinity();
                    }
                }else{
                    b.disorientation = std::numeric_limits<double>::infinity();
                }
            }
            progress.advance(r.size());
        }, tbb::auto_partitioner{});

        tbb::parallel_sort(_neighborBonds.begin(), _neighborBonds.end(), [](const NeighborBond& x, const NeighborBond& y){
            if(x.disorientation != y.disorientation) return x.disorientation < y.disorientation;
            if(x.a != y.a) return x.a < y.a;
            return x.b < y.b;
        });
    }

    double calculate_disorientation(int structureType, Quaternion& qa, const Quaternion& qb){
        double qa_norm = qa.norm();
        double qb_norm = qb.norm();
        double qtarget[4] = { qa.w()/qa_norm, qa.x()/qa_norm, qa.y()/qa_norm, qa.z()/qa_norm };
        double q[4] = { qb.w()/qb_norm, qb.x()/qb_norm, qb.y()/qb_norm, qb.z()/qb_norm };

        int type = 0;
        if(structureType == StructureType::OTHER){
            return std::numeric_limits<double>::max();
        }else if(structureType == StructureType::FCC){
            type = PTM_MATCH_FCC;
        }else if(structureType == StructureType::HCP){
            type = PTM_MATCH_HCP;
        }else if(structureType == StructureType::BCC){
            type = PTM_MATCH_BCC;
        }else if(structureType == StructureType::SC){
            type = PTM_MATCH_SC;
        }else if(structureType == StructureType::CUBIC_DIAMOND){
            type = PTM_MATCH_DCUB;
        }else if(structureType == StructureType::HEX_DIAMOND){
            type = PTM_MATCH_DHEX;
        }else if(structureType == StructureType::GRAPHENE){
            type = PTM_MATCH_GRAPHENE;
        }

        double disorientation = (double) ptm_map_and_calculate_disorientation(type, qtarget, q);

        qa.w() += q[0] * qb_norm;
        qa.x() += q[1] * qb_norm;
        qa.y() += q[2] * qb_norm;
        qa.z() += q[3] * qb_norm;

        return disorientation;
    }

    void node_pair_sampling_clustering(Graph& graph, std::vector<Quaternion>& qsum){
        double totalWeight = 1;

        size_t progressVal = 0;
        // Every contraction or isolation retires one node, so retired nodes measure the progress.
        const size_t initialNodes = graph.num_nodes();
        ProgressReporter progress("  Graph clustering", initialNodes);
        std::vector<size_t> chain;
        while(graph.num_nodes()){
            size_t node = graph.next_node();

            chain.push_back(node);
            while(!chain.empty()){
                progress.update(initialNodes - graph.num_nodes());

                size_t a = chain.back();
                chain.pop_back();

                auto [d, b] = graph.nearestNeighbor(a);
                if(b == std::numeric_limits<size_t>::max()){
                    graph.remove_node(a);
                }else if(!chain.empty()){
                    size_t c = chain.back();
                    chain.pop_back();

                    if(b == c){
                        size_t parent = graph.contract_edge(a, b);
                        size_t child = (parent == a) ? b : a;

                        double disorientation = calculate_disorientation(_adjustedStructureTypes[parent], qsum[parent], qsum[child]);
                        _dendrogram.emplace_back(parent, child, d / totalWeight, disorientation, 1, qsum[parent]);
                    }else{
                        chain.push_back(c);
                        chain.push_back(a);
                        chain.push_back(b);
                    }
                }else{
                    chain.push_back(a);
                    chain.push_back(b);
                }
            }
        }
    }

    void minimum_spanning_tree_clustering(std::vector<Quaternion>& qsum, DisjointSet& uf){
        for(const NeighborBond& edge : _neighborBonds){
            if(edge.disorientation < _misorientationThresholdDeg){
                size_t pa = uf.find(edge.a);
                size_t pb = uf.find(edge.b);
                if(pa != pb && isCrystallineBond(edge)){
                    size_t parent = uf.merge(pa, pb);
                    size_t child = (parent == pa) ? pb : pa;
                    double disorientation = calculate_disorientation(_adjustedStructureTypes[parent], qsum[parent], qsum[child]);
                    _dendrogram.emplace_back(parent, child, edge.disorientation, disorientation, 1, qsum[parent]);
                }
            }
        }
    }

    void determineMergeSequence(){
        Graph graph(_numParticles, _neighborBonds.size());
        if(_algorithmType == MergeAlgorithm::GraphClusteringAutomatic || _algorithmType == MergeAlgorithm::GraphClusteringManual){
            for(auto edge : _neighborBonds){
                if(joinsMergeGraph(edge)){
                    double weight = calculateGraphWeight(edge.disorientation);
                    graph.add_edge(edge.a, edge.b, weight);
                }
            }
        }

        std::vector<Quaternion> qsum(_adjustedOrientations.cbegin(), _adjustedOrientations.cend());
        DisjointSet uf(_numParticles);
        _dendrogram.resize(0);

        if(_algorithmType == MergeAlgorithm::GraphClusteringAutomatic || _algorithmType == MergeAlgorithm::GraphClusteringManual){
            node_pair_sampling_clustering(graph, qsum);
        }else{
            minimum_spanning_tree_clustering(qsum, uf);
        }

        std::sort(_dendrogram.begin(), _dendrogram.end(), [](const DendrogramNode& a, const DendrogramNode& b){
            if(a.distance != b.distance) return a.distance < b.distance;
            if(a.a != b.a) return a.a < b.a;
            return a.b < b.b;
        });

        size_t numPlot = 0;
        uf.clear();
        for(DendrogramNode& node : _dendrogram){
            size_t sa = uf.nodesize(uf.find(node.a));
            size_t sb = uf.nodesize(uf.find(node.b));
            size_t dsize = std::min(sa, sb);

            node.merge_size = 2. / (1. / sa + 1. / sb);
            uf.merge(node.a, node.b);

            node.size = dsize;
            if(dsize >= _minPlotSize){
                numPlot++;
            }
        }

        if(_algorithmType == MergeAlgorithm::GraphClusteringAutomatic || _algorithmType == MergeAlgorithm::GraphClusteringManual){
            std::vector<double> mergeDistanceArray;
            std::vector<double> mergeSizeArray;

            mergeDistanceArray.reserve(numPlot);
            mergeSizeArray.reserve(numPlot);

            for(const DendrogramNode& node : _dendrogram){
                if(node.size >= _minPlotSize){
                    mergeDistanceArray.push_back(std::log(node.distance));
                    mergeSizeArray.push_back(node.size);
                }
            }

            auto regressor = Regressor(_dendrogram);
            _suggestedMergingThreshold = regressor.calculate_threshold(_dendrogram, 1.5);

            numPlot = 0;
            for(auto y : regressor.ys){
                numPlot += (y > 0) ? 1 : 0;
            }

            std::vector<double> logMergeSizeArray;
            std::vector<double> logMergeDistanceArray;

            logMergeSizeArray.reserve(numPlot);
            logMergeDistanceArray.reserve(numPlot);

            for(size_t i = 0; i < regressor.residuals.size(); i++){
                if(regressor.ys[i] > 0){
                    logMergeSizeArray.push_back(regressor.xs[i]);
                    logMergeDistanceArray.push_back(regressor.ys[i]);
                }
            }
        }
    }

private:
    static constexpr double _misorientationThresholdDeg = 4.0;
    const size_t _minPlotSize = 20;

    MergeAlgorithm _algorithmType;
    bool _handleBoundaries;
    size_t _numParticles;

    std::shared_ptr<ParticleProperty> _positions;
    std::shared_ptr<ParticleProperty> _structuresProperty;
    std::shared_ptr<ParticleProperty> _orientationsProperty;
    std::shared_ptr<ParticleProperty> _correspondencesProperty;

    const SimulationCell _simCell;
    bool _outputBonds;

    std::vector<NeighborBond> _neighborBonds;
    std::vector<StructureType> _adjustedStructureTypes;
    std::vector<Quaternion> _adjustedOrientations;

    std::vector<DendrogramNode> _dendrogram;
    double _suggestedMergingThreshold = 0.0;
};

class GrainSegmentationEngine2{
public:
    struct GrainInfo{
        int id;
        size_t size;
        Quaternion orientation;
        int structureType;
    };

    GrainSegmentationEngine2(
        std::shared_ptr<const GrainSegmentationEngine1> engine1,
        double mergingThreshold,
        bool adoptOrphanAtoms,
        size_t minGrainAtomCount,
        bool colorParticlesByGrain
    )
    : _engine1(std::move(engine1))
    , _numParticles(_engine1 ? _engine1->structuresProperty()->size() : 0)
    , _mergingThreshold(mergingThreshold)
    , _adoptOrphanAtoms(adoptOrphanAtoms)
    , _minGrainAtomCount(minGrainAtomCount)
    , _colorParticlesByGrain(colorParticlesByGrain)
    {
        _atomClusters = std::make_shared<ParticleProperty>(_numParticles, DataType::Int, 1, 0, false);
    }

    // The merges up to the merging threshold of the dendrogram applied, out of all its merges; the threshold as
    // the line names it (in degrees for the minimum spanning tree).
    struct MergesApplied{
        size_t applied;
        size_t merges;
        MergeAlgorithm algorithm;
        double threshold;
    };

    // The atoms of no grain that joined a neighboring one, out of the atoms of no grain.
    struct OrphanAdoption{
        size_t adoptedAtoms;
        size_t orphanAtoms;
    };

    void perform(GrainSegmentationLog& log);

    size_t grainCount() const{
        return _grainCount;
    }

    const std::vector<GrainInfo>& grains() const{
        return _grains;
    }

    std::shared_ptr<ParticleProperty> atomClusters() const{
        return _atomClusters;
    }

    size_t atomCount() const{
        return _numParticles;
    }

    size_t minGrainAtomCount() const{
        return _minGrainAtomCount;
    }

    bool adoptsOrphanAtoms() const{
        return _adoptOrphanAtoms;
    }

private:
    struct PQNode{
        int cluster;
        size_t particleIndex;
        double length;
    };

    struct PQCompareLength{
        bool operator()(const PQNode& a, const PQNode& b) const{
            return a.length > b.length;
        }
    };

    OrphanAdoption mergeOrphanAtoms(){
        using NeighborBond = GrainSegmentationEngine1::NeighborBond;
        size_t orphanAtoms = _numParticles;
        for(const auto& grain : _grains){
            orphanAtoms -= grain.size;
        }
        size_t adoptedAtoms = 0;

        std::vector<NeighborBond> noncrystallineBonds;
        for(auto nb : _engine1->neighborBonds()){
            if(_atomClusters->getInt((size_t) nb.a) == 0 || _atomClusters->getInt((size_t) nb.b) == 0){
                noncrystallineBonds.push_back(nb);

                std::swap(nb.a, nb.b);
                noncrystallineBonds.push_back(nb);
            }
        }

        std::sort(noncrystallineBonds.begin(), noncrystallineBonds.end(),
            [](const NeighborBond& a, const NeighborBond& b){ return a.a < b.a; });

        std::priority_queue<PQNode, std::vector<PQNode>, PQCompareLength> pq;

        for(const auto& bond : _engine1->neighborBonds()){
            int clusterA = _atomClusters->getInt((size_t) bond.a);
            int clusterB = _atomClusters->getInt((size_t) bond.b);

            if(clusterA != 0 && clusterB == 0){
                pq.push({ clusterA, bond.b, bond.length });
            }else if(clusterA == 0 && clusterB != 0){
                pq.push({ clusterB, bond.a, bond.length });
            }
        }

        ProgressReporter progress("  Adopting orphan atoms", orphanAtoms);
        while(!pq.empty()){
            PQNode node = pq.top();
            pq.pop();

            if(_atomClusters->getInt(node.particleIndex) != 0)
                continue;

            _atomClusters->setInt(node.particleIndex, node.cluster);
            progress.update(++adoptedAtoms);

            if(node.cluster >= 1 && (size_t) node.cluster <= _grains.size())
                _grains[node.cluster - 1].size++;

            NeighborBond key{ node.particleIndex, 0, 0.0, 0.0 };
            auto bondsRange = std::equal_range(noncrystallineBonds.begin(), noncrystallineBonds.end(), key,
                [](const NeighborBond& a, const NeighborBond& b){ return a.a < b.a; });

            for(auto it = bondsRange.first; it != bondsRange.second; ++it){
                const NeighborBond& bond = *it;

                size_t neighborIndex = bond.b;
                if(neighborIndex == std::numeric_limits<size_t>::max()) break;
                if(_atomClusters->getInt(neighborIndex) != 0) continue;

                pq.push({ node.cluster, neighborIndex, node.length + bond.length });
            }
        }
        return {adoptedAtoms, orphanAtoms};
    }

    std::shared_ptr<const GrainSegmentationEngine1> _engine1;
    size_t _numParticles = 0;

    double _mergingThreshold = 0.0;
    bool _adoptOrphanAtoms = false;
    size_t _minGrainAtomCount = 1;
    bool _colorParticlesByGrain = false;
    size_t _grainCount = 0;

    std::vector<GrainInfo> _grains;
    std::shared_ptr<ParticleProperty> _atomClusters;
};

}
