#include <volt/grain_segmentation_log.h>

#include <volt/analysis/structure_analysis.h>
#include <volt/core/type_counts.h>
#include <volt/structures/crystal_structure_types.h>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <array>
#include <string>
#include <vector>

namespace Volt{

void GrainSegmentationLog::dump(const StructureAnalysis& identified){
    const StructureContext& context = identified.context();
    const std::size_t atoms = context.atomCount();
    constexpr int K = static_cast<int>(StructureType::NUM_STRUCTURE_TYPES);
    std::array<std::size_t, K> counts{};
    for(std::size_t i = 0; i < atoms; ++i){
        const int type = context.structureTypes->getInt(i);
        ++counts[(0 <= type && type < K) ? type : 0];
    }
    TypeCounts line{"PTM structure identification done in " + phaseElapsed(), {}, atoms};
    for(const StructureType type : {StructureType::FCC, StructureType::HCP, StructureType::BCC, StructureType::ICO, StructureType::SC}){
        line.counts.emplace_back(structureTypeName(type), counts[type]);
    }
    line.counts.emplace_back("Other", counts[StructureType::OTHER]);
    dump(line);
    _engine.restart();
}

void GrainSegmentationLog::dump(const GrainSegmentationEngine1::NeighborBonds& neighborBonds){
    spdlog::info("  Neighbor bonds: {} bonds between {} atoms in {}",
        neighborBonds.bonds.size(), neighborBonds.atoms, phaseElapsed());
}

void GrainSegmentationLog::dump(const GrainSegmentationEngine1::CoherentInterfaces& interfaces){
    spdlog::info("  Coherent interfaces: rotated {} atoms onto their parent phase (FCC/HCP parent {}) in {}",
        interfaces.rotatedAtoms, structureTypeName(interfaces.fccParent), phaseElapsed());
}

void GrainSegmentationLog::dump(const GrainSegmentationEngine1::BondDisorientations& disorientations){
    spdlog::info("  Disorientations: {} bonds computed and sorted in {}", disorientations.bonds.size(), phaseElapsed());
}

std::pair<std::size_t, std::size_t> GrainSegmentationLog::mergeGraphSize(const GrainSegmentationEngine1& engine){
    std::vector<bool> inGraph(engine.structuresProperty()->size(), false);
    std::size_t atoms = 0;
    std::size_t bonds = 0;
    for(const auto& bond : engine.neighborBonds()){
        if(!engine.joinsMergeGraph(bond)) continue;
        ++bonds;
        for(const std::size_t atom : {bond.a, bond.b}){
            if(!inGraph[atom]){
                inGraph[atom] = true;
                ++atoms;
            }
        }
    }
    return {atoms, bonds};
}

void GrainSegmentationLog::dump(const GrainSegmentationEngine1& engine){
    const std::string mergeElapsed = phaseElapsed();
    const std::string engineElapsed = _engine.elapsed();
    const double threshold = GrainSegmentationEngine1::misorientationThresholdDeg();
    if(engine.algorithmType() == MergeAlgorithm::MinimumSpanningTree){
        spdlog::info("  Merge sequence: {} merges from the minimum spanning tree of bonds below {} degrees in {}",
            engine.dendrogram().size(), threshold, mergeElapsed);
        spdlog::info("GrainSegmentationEngine1 complete in {}.", engineElapsed);
    }else{
        const auto [graphAtoms, graphBonds] = mergeGraphSize(engine);
        spdlog::info("  Merge sequence: {} merges from graph clustering of {} atoms over {} bonds below {} degrees in {}",
            engine.dendrogram().size(), graphAtoms, graphBonds, threshold, mergeElapsed);
        spdlog::info("GrainSegmentationEngine1 complete in {}. Suggested merging threshold: {:.4f}",
            engineElapsed, engine.suggestedMergingThreshold());
    }
    _engine.restart();
}

void GrainSegmentationLog::dump(const GrainSegmentationEngine2::MergesApplied& merges){
    if(merges.algorithm == MergeAlgorithm::MinimumSpanningTree){
        spdlog::info("  Applied {} of {} merges with disorientation up to {} degrees",
            merges.applied, merges.merges, merges.threshold);
    }else{
        spdlog::info("  Applied {} of {} merges up to merging threshold {:.4f} ({})",
            merges.applied, merges.merges, merges.threshold,
            merges.algorithm == MergeAlgorithm::GraphClusteringAutomatic ? "suggested" : "manual");
    }
    // The line names no time, and the next one times the orphan adoption alone.
    phaseMilliseconds();
}

void GrainSegmentationLog::dump(const GrainSegmentationEngine2::OrphanAdoption& adoption){
    spdlog::info("  Orphan adoption: {} of {} unassigned atoms joined a neighboring grain in {}",
        adoption.adoptedAtoms, adoption.orphanAtoms, phaseElapsed());
}

void GrainSegmentationLog::dump(const GrainSegmentationEngine2& engine){
    const std::string elapsed = _engine.elapsed();
    const std::size_t atoms = engine.atomCount();
    // Grains are few, so a copy of their sizes for the median is cheap.
    std::vector<std::size_t> grainSizes;
    grainSizes.reserve(engine.grains().size());
    std::size_t groupedAtoms = 0;
    for(const auto& grain : engine.grains()){
        grainSizes.push_back(grain.size);
        groupedAtoms += grain.size;
    }
    const std::size_t ungroupedAtoms = atoms > groupedAtoms ? atoms - groupedAtoms : 0;
    if(grainSizes.empty()){
        spdlog::warn("Found no grains with at least {} atoms ({}); all {} atoms are unassigned",
            engine.minGrainAtomCount(), elapsed, atoms);
        return;
    }
    const auto [smallest, largest] = std::minmax_element(grainSizes.begin(), grainSizes.end());
    const std::size_t smallestSize = *smallest;
    const std::size_t largestSize = *largest;
    std::nth_element(grainSizes.begin(), grainSizes.begin() + grainSizes.size() / 2, grainSizes.end());
    spdlog::info("Found {} grains in {}: largest {} atoms, median {}, smallest {}; {} atoms ({:.1f}%) unassigned",
        engine.grainCount(), elapsed, largestSize, grainSizes[grainSizes.size() / 2], smallestSize,
        ungroupedAtoms, percentOf(ungroupedAtoms, atoms));
    if(engine.adoptsOrphanAtoms() && ungroupedAtoms > 0){
        spdlog::warn("{} atoms have no bonded path to any grain and stay unassigned", ungroupedAtoms);
    }
}

}
