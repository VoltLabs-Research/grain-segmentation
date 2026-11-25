#pragma once

#include <volt/core/run_log.h>
#include <volt/grain_segmentation_engine.h>

#include <cstddef>
#include <utility>

namespace Volt{

class StructureAnalysis;

// The execution log of one grain segmentation. Each phase hands it what it built, and the log writes what that holds
// and how long the phase took. Scripts parse these lines, so their wording stays as it is.
class GrainSegmentationLog : public RunLog{
public:
    using RunLog::dump;

    // The structures PTM identified, counted by type.
    void dump(const StructureAnalysis& identified);

    void dump(const GrainSegmentationEngine1::NeighborBonds& neighborBonds);
    void dump(const GrainSegmentationEngine1::CoherentInterfaces& interfaces);
    void dump(const GrainSegmentationEngine1::BondDisorientations& disorientations);
    // The merge sequence, the engine's last phase, then the engine as a whole.
    void dump(const GrainSegmentationEngine1& engine);

    void dump(const GrainSegmentationEngine2::MergesApplied& merges);
    void dump(const GrainSegmentationEngine2::OrphanAdoption& adoption);
    // The grains, timed over the whole engine.
    void dump(const GrainSegmentationEngine2& engine);

private:
    // The atoms and the bonds of the graph the merge sequence clustered.
    static std::pair<std::size_t, std::size_t> mergeGraphSize(const GrainSegmentationEngine1& engine);

    // Since the engine whose closing line comes next started: the dump before each engine restarts it.
    Stopwatch _engine;
};

}
