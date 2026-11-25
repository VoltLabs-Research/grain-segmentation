#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>
#include <nlohmann/json.hpp>
#include <volt/core/particle_property.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/ptm_structure_analysis.h>
#include <volt/grain_segmentation_engine.h>
#include <ptm_constants.h>
#include <string>
#include <vector>

namespace Volt{
using json = nlohmann::json;

class GrainSegmentationLog;

class GrainSegmentationService{
public:
    GrainSegmentationService();

    void setRMSD(float rmsd);
    void setMinGrainAtomCount(int minGrainAtomCount);
    void setAdoptOrphanAtoms(bool adoptOrphanAtoms);
    void setHandleCoherentInterfaces(bool handleCoherentInterfaces);
    void setMergeAlgorithm(MergeAlgorithm mergeAlgorithm);
    void setMergingThreshold(double mergingThreshold);
    void setOutputBonds(bool outputBonds);

    json compute(
        const LammpsParser::Frame &frame,
        const std::string &outputFilename = ""
    );

private:
    static constexpr int kSearchedPtmStructures =
        PTM_CHECK_SC | PTM_CHECK_FCC | PTM_CHECK_HCP | PTM_CHECK_ICO | PTM_CHECK_BCC;

    float _rmsd;

    bool _adoptOrphanAtoms;
    int _minGrainAtomCount;
    bool _handleCoherentInterfaces;
    bool _outputBonds;
    MergeAlgorithm _mergeAlgorithm;
    double _mergingThreshold;

    json performGrainSegmentation(
        const LammpsParser::Frame &frame,
        const std::vector<int>& structureTypes,
        const std::vector<PtmLocalAtomState>& ptmStates,
        const std::string& outputFile,
        GrainSegmentationLog& log
    );

    // The grain table: each grain with its size, structure, orientation and center of mass.
    static json grainTable(
        const LammpsParser::Frame& frame,
        const GrainSegmentationEngine1& engine1,
        const GrainSegmentationEngine2& engine2,
        const std::vector<int>& grainIds
    );

    // The atoms table but its per-atom properties: the atoms grouped by structure, with their grains.
    static json atomsByStructure(
        const LammpsParser::Frame& frame,
        const std::vector<int>& structureTypes,
        const std::vector<PtmLocalAtomState>& ptmStates,
        const std::vector<int>& grainIds
    );

    static json perAtomProperties(
        const LammpsParser::Frame& frame,
        const std::vector<int>& structureTypes,
        const std::vector<PtmLocalAtomState>& ptmStates,
        const std::vector<int>& grainIds
    );
};

}
