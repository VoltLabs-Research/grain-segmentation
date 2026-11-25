#include <volt/grain_segmentation_engine.h>
#include <volt/grain_segmentation_log.h>

#include <algorithm>
#include <cmath>
#include <numeric>

namespace Volt{

void GrainSegmentationEngine1::perform(GrainSegmentationLog& log){
    // The graph clustering of the merge sequence runs serially over every atom; the minimum
    // spanning tree is cheaper.
    ProgressSteps steps({
        {"Building neighbor bonds", 1.0},
        {"Rotating coherent interface atoms", 1.0, _handleBoundaries},
        {"Computing bond disorientations", 1.0},
        {"Determining the grain merge sequence", _algorithmType == MergeAlgorithm::MinimumSpanningTree ? 1.0 : 5.0},
    });
    createNeighborBonds();
    log.dump(NeighborBonds{_neighborBonds, _numParticles});

    steps.next();
    if(_handleBoundaries){
        log.dump(rotateInterfaceAtoms());
    }

    steps.next();
    computeDisorientationAngles();
    log.dump(BondDisorientations{_neighborBonds});

    steps.next();
    determineMergeSequence();

    _positions.reset();
}

void GrainSegmentationEngine2::perform(GrainSegmentationLog& log){
    if(!_engine1) return;

    const auto& dendro = _engine1->dendrogram();

    double thr = _mergingThreshold;
    if(_engine1->algorithmType() == MergeAlgorithm::GraphClusteringAutomatic){
        thr = _engine1->suggestedMergingThreshold();
    }
    if(_engine1->algorithmType() == MergeAlgorithm::MinimumSpanningTree){
        thr = std::log(thr);
    }

    DisjointSet uf(_numParticles);
    std::vector<Quaternion> meanQ(_engine1->orientationsProperty()->size());
    const double* qptr = _engine1->orientationsProperty()->dataDouble();

    for(size_t i = 0; i < _numParticles; ++i){
        meanQ[i] = Quaternion(qptr[4 * i + 0], qptr[4 * i + 1], qptr[4 * i + 2], qptr[4 * i + 3]);
    }

    // The dendrogram is sorted by distance, so the merges up to the threshold are a prefix of it.
    const auto applied = std::ranges::find_if(dendro, [thr](const GrainSegmentationEngine1::DendrogramNode& node){
        return std::log(node.distance) > thr;
    });
    for(auto node = dendro.begin(); node != applied; ++node){
        uf.merge(node->a, node->b);
        size_t p = uf.find(node->a);
        meanQ[p] = node->orientation;
    }
    log.dump(MergesApplied{
        static_cast<size_t>(applied - dendro.begin()),
        dendro.size(),
        _engine1->algorithmType(),
        _engine1->algorithmType() == MergeAlgorithm::MinimumSpanningTree ? _mergingThreshold : thr
    });

    std::vector<size_t> rep2id(_numParticles, 0);
    size_t nextId = 1;
    for(size_t i = 0; i < _numParticles; ++i){
        if(uf.find(i) == i){
            if(uf.nodesize(i) < _minGrainAtomCount || _engine1->structuresProperty()->getInt(i) == (int)StructureType::OTHER){
                rep2id[i] = 0;
            }else{
                rep2id[i] = nextId++;
            }
        }
    }

    for(size_t i = 0; i < _numParticles; ++i){
        size_t rep = uf.find(i);
        _atomClusters->setInt(i, (int) rep2id[rep]);
    }

    _grainCount = nextId - 1;
    _grains.clear();
    _grains.reserve(_grainCount);

    for(size_t rep = 0; rep < _numParticles; ++rep){
        if(uf.find(rep) == rep){
            int gid = (int) rep2id[rep];
            if(gid > 0){
                Quaternion q = meanQ[rep].normalized();
                int st = _engine1->structuresProperty()->getInt(rep);
                _grains.emplace_back(GrainInfo{gid, uf.nodesize(rep), q, st});
            }
        }
    }

    if(_grainCount >= 1){
        std::vector<size_t> mapping(_grainCount);
        std::iota(mapping.begin(), mapping.end(), size_t(0));
        std::sort(mapping.begin(), mapping.end(), [&](size_t a, size_t b){
            return _grains[a].size > _grains[b].size;
        });

        std::vector<GrainInfo> reordered;
        reordered.reserve(_grainCount);
        for(size_t p = 0; p < _grainCount; ++p){
            GrainInfo g = _grains[mapping[p]];
            g.id = (int)(p + 1);
            reordered.push_back(g);
        }
        _grains = std::move(reordered);

        std::vector<size_t> inverseMapping(_grainCount + 1);
        inverseMapping[0] = 0;
        for(size_t i = 1; i <= _grainCount; ++i){
            inverseMapping[mapping[i - 1] + 1] = i;
        }

        for(size_t i = 0; i < _numParticles; ++i){
            size_t oldId = (size_t) _atomClusters->getInt(i);
            _atomClusters->setInt(i, (int) inverseMapping[oldId]);
        }

        if(_adoptOrphanAtoms)
            log.dump(mergeOrphanAtoms());
    }
}

}
