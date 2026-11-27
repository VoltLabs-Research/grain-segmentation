#include <volt/grain_segmentation_service.h>
#include <volt/grain_segmentation_log.h>
#include <volt/analysis/ptm_structure_analysis.h>
#include <volt/analysis/structure_identification_export.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/core/progress.h>
#include <volt/plugin/output_files.h>
#include <volt/utilities/json_utils.h>
#include <volt/utilities/written_table.h>
#include <map>
#include <cstdint>
#include <utility>

namespace Volt{

using namespace Volt::Particles;

GrainSegmentationService::GrainSegmentationService()
    : _rmsd(0.10f),
      _adoptOrphanAtoms(true),
      _minGrainAtomCount(100),
      _handleCoherentInterfaces(true),
      _outputBonds(false),
      _mergeAlgorithm(MergeAlgorithm::GraphClusteringAutomatic),
      _mergingThreshold(0.0){}

void GrainSegmentationService::setRMSD(float rmsd){
    _rmsd = rmsd;
}

void GrainSegmentationService::setMinGrainAtomCount(int minGrainAtomCount){
    _minGrainAtomCount = minGrainAtomCount;
}

void GrainSegmentationService::setAdoptOrphanAtoms(bool adoptOrphanAtoms){
    _adoptOrphanAtoms = adoptOrphanAtoms;
}

void GrainSegmentationService::setHandleCoherentInterfaces(bool handleCoherentInterfaces){
    _handleCoherentInterfaces = handleCoherentInterfaces;
}

void GrainSegmentationService::setMergeAlgorithm(MergeAlgorithm mergeAlgorithm){
    _mergeAlgorithm = mergeAlgorithm;
}

void GrainSegmentationService::setMergingThreshold(double mergingThreshold){
    _mergingThreshold = mergingThreshold;
}

void GrainSegmentationService::setOutputBonds(bool outputBonds){
    _outputBonds = outputBonds;
}

json GrainSegmentationService::compute(const LammpsParser::Frame &frame, const std::string &outputFilename){
    if(frame.simulationCell.is2D())
        return AnalysisResult::failure("The grain segmentation modifier does not support 2d simulation cells.");

    GrainSegmentationLog log;
    FrameAdapter::PreparedAnalysisInput prepared;
    std::string frameError;
    if(!FrameAdapter::prepareAnalysisInput(frame, prepared, &frameError))
        return AnalysisResult::failure(frameError);

    auto positions = std::move(prepared.positions);

    std::vector<Matrix3> preferredOrientations;
    preferredOrientations.push_back(Matrix3::Identity());

    auto structuretypes = std::make_unique<ParticleProperty>(frame.natoms, DataType::Int, 1, 0, true);
    AnalysisContext context(
        positions.get(),
        frame.simulationCell,
        LATTICE_OTHER,
        structuretypes.get(),
        std::move(preferredOrientations)
    );

    auto structureAnalysis = std::make_unique<StructureAnalysis>(context);
    auto ptmStates = std::make_shared<std::vector<PtmLocalAtomState>>();
    // The segmentation's outputs build a JSON object per atom.
    ProgressSteps steps({
        {"Identifying local structures with PTM", 1.0},
        {"Segmenting grains", 8.0},
    });
    determineLocalStructuresWithPTM(*structureAnalysis, _rmsd, ptmStates, nullptr, kSearchedPtmStructures);
    computeMaximumNeighborDistanceFromPTM(*structureAnalysis);
    log.dump(*structureAnalysis);

    std::vector<int> extractedStructureTypes;
    extractedStructureTypes.reserve(frame.natoms);
    for(int i = 0; i < frame.natoms; i++){
        extractedStructureTypes.push_back(structureAnalysis->context().structureTypes->getInt(i));
    }
    steps.next();

    if(!outputFilename.empty()){
        return performGrainSegmentation(frame, extractedStructureTypes, *ptmStates, outputFilename, log);
    }

    return AnalysisResult::failure("No output filename specified");
}

json GrainSegmentationService::performGrainSegmentation(
    const LammpsParser::Frame &frame,
    const std::vector<int>& structureTypes,
    const std::vector<PtmLocalAtomState>& ptmStates,
    const std::string &outputFile,
    GrainSegmentationLog& log
){
    try{
        if(ptmStates.size() < static_cast<size_t>(frame.natoms)){
            return AnalysisResult::failure("Grain segmentation requires PTM orientation state for all atoms.");
        }

        // The atoms table serializes both JSON listings.
        ProgressSteps steps({
            {"Running grain segmentation engine 1", 7.0},
            {"Running grain segmentation engine 2", 1.0},
            {"Building the grain table", 3.0},
            {"Collecting per-atom properties", 3.0},
            {"Writing the atoms table", 2.0},
        });

        auto positions = std::make_shared<ParticleProperty>(frame.natoms, ParticleProperty::PositionProperty, 0, true);
        for(size_t i = 0; i < frame.positions.size() && i < static_cast<size_t>(frame.natoms); i++){
            positions->setPoint3(i, frame.positions[i]);
        }

        auto structures = std::make_shared<ParticleProperty>(frame.natoms, DataType::Int, 1, 0, false);
        for(size_t i = 0; i < structureTypes.size(); i++){
            structures->setInt(i, structureTypes[i]);
        }

        auto orientations = std::make_shared<ParticleProperty>(frame.natoms, DataType::Double, 4, 0, false);
        for(size_t i = 0; i < static_cast<size_t>(frame.natoms); i++){
            const Quaternion q = ptmStates[i].orientation.normalized();
            orientations->setDoubleComponent(i, 0, q.x());
            orientations->setDoubleComponent(i, 1, q.y());
            orientations->setDoubleComponent(i, 2, q.z());
            orientations->setDoubleComponent(i, 3, q.w());
        }

        auto correspondences = std::make_shared<ParticleProperty>(frame.natoms, DataType::Int64, 1, 0, false);
        for(size_t i = 0; i < static_cast<size_t>(frame.natoms); ++i){
            correspondences->setInt64(i, 0);
        }

        auto engine1 = std::make_shared<GrainSegmentationEngine1>(
            positions,
            structures,
            orientations,
            correspondences,
            &frame.simulationCell,
            _mergeAlgorithm,
            _handleCoherentInterfaces,
            _outputBonds
        );
        engine1->perform(log);
        log.dump(*engine1);

        steps.next();
        GrainSegmentationEngine2 engine2(
            engine1,
            _mergingThreshold,
            _adoptOrphanAtoms,
            static_cast<size_t>(_minGrainAtomCount),
            true
        );
        engine2.perform(log);
        log.dump(engine2);

        steps.next();
        auto atomClusters = engine2.atomClusters();
        std::vector<int> grainIds(frame.natoms, 0);
        for(size_t i = 0; i < static_cast<size_t>(frame.natoms); i++){
            grainIds[i] = atomClusters->getInt(i);
        }
        json result = grainTable(frame, *engine1, engine2, grainIds);
        log.dump(JsonUtils::writeSummaryTable(result, outputFile, "grains", "Grains"));
        json atomsTable = atomsByStructure(frame, structureTypes, ptmStates, grainIds);

        steps.next();
        atomsTable["per-atom-properties"] = perAtomProperties(frame, structureTypes, ptmStates, grainIds);

        steps.next();
        const AtomsTable atoms{outputPath(outputFile, "atoms")};
        JsonUtils::writeJsonToParquet(atomsTable, atoms.path, false);
        log.dump(atoms);

        return result;
    }catch(const std::exception& e){
        return AnalysisResult::failure(std::string("Grain segmentation failed: ") + e.what());
    }
}

json GrainSegmentationService::grainTable(
    const LammpsParser::Frame& frame,
    const GrainSegmentationEngine1& engine1,
    const GrainSegmentationEngine2& engine2,
    const std::vector<int>& grainIds
){
    const SimulationCell& cell = frame.simulationCell;
    std::map<int, Point3> grainReference;
    std::map<int, Vector3> grainOffsetSum;
    std::map<int, int> grainAtomCount;
    for(int i = 0; i < frame.natoms; i++){
        int gid = grainIds[i];
        if(gid == 0) continue;
        if(i >= static_cast<int>(frame.positions.size())) continue;
        const Point3& p = frame.positions[i];
        auto it = grainReference.find(gid);
        if(it == grainReference.end()){
            grainReference[gid] = p;
            grainOffsetSum[gid] = Vector3(Vector3::Zero{});
            grainAtomCount[gid] = 1;
        }else{
            grainOffsetSum[gid] += cell.unwrapVector(p - it->second);
            grainAtomCount[gid]++;
        }
    }

    json grainsArray = json::array();
    for(const auto &grain : engine2.grains()){
        json grainInfo;
        grainInfo["id"] = grain.id;
        grainInfo["size"] = grain.size;
        grainInfo["structure_type"] = grain.structureType;
        grainInfo["structure_name"] = StructureIdentificationExport::structureTypeNameForExport(grain.structureType);
        grainInfo["orientation"] = grain.orientation;
        if(grainAtomCount.count(grain.id) && grainAtomCount[grain.id] > 0){
            int cnt = grainAtomCount[grain.id];
            grainInfo["pos"] = cell.wrapPoint(grainReference[grain.id] + (grainOffsetSum[grain.id] / cnt));
        } else {
            grainInfo["pos"] = {0.0, 0.0, 0.0};
        }
        grainsArray.push_back(grainInfo);
    }

    json result;
    result["main_listing"] = {
        { "total_grains", static_cast<int>(engine2.grainCount()) },
        { "merging_threshold", engine1.suggestedMergingThreshold() }
    };
    result["sub_listings"] = { { "grains", grainsArray } };
    return result;
}

json GrainSegmentationService::atomsByStructure(
    const LammpsParser::Frame& frame,
    const std::vector<int>& structureTypes,
    const std::vector<PtmLocalAtomState>& ptmStates,
    const std::vector<int>& grainIds
){
    const StructureIdentificationExport::AtomsByStructure structures = StructureIdentificationExport::groupAtomsByStructure(
        structureTypes, [&](std::size_t atomIndex){ return grainIds[atomIndex] != 0; });

    ProgressReporter structureProgress("  Serializing atoms by structure", static_cast<std::uint64_t>(frame.natoms));
    size_t serializedAtoms = 0;
    return StructureIdentificationExport::atomsByStructureTable(frame, structures, [&](json& atom, std::size_t atomIndex){
        structureProgress.update(++serializedAtoms);
        const PtmLocalAtomState& state = ptmStates[atomIndex];
        atom["grain_id"] = grainIds[atomIndex];
        atom["orientation"] = state.orientation.normalized();
        atom["ptm_valid"] = state.valid;
        if(state.valid){
            atom["rmsd"] = state.rmsd;
        }
    });
}

json GrainSegmentationService::perAtomProperties(
    const LammpsParser::Frame& frame,
    const std::vector<int>& structureTypes,
    const std::vector<PtmLocalAtomState>& ptmStates,
    const std::vector<int>& grainIds
){
    constexpr int K = static_cast<int>(StructureType::NUM_STRUCTURE_TYPES);
    json perAtom = json::array();
    ProgressReporter perAtomProgress("  Serializing per-atom properties", static_cast<std::uint64_t>(frame.natoms));
    for(size_t i = 0; i < static_cast<size_t>(frame.natoms); ++i){
        perAtomProgress.update(i + 1);
        const int raw = structureTypes[i];
        const int st = (0 <= raw && raw < K) ? raw : 0;
        const PtmLocalAtomState& state = ptmStates[i];
        const int grainId = grainIds[i];
        json p = {
            {"id", frame.ids[i]},
            {"structure_id", st},
            {"structure_name", StructureIdentificationExport::structureTypeNameForExport(st)},
            {"grain_id", grainId},
            {"orientation", state.orientation.normalized()},
            {"ptm_valid", state.valid}
        };
        if(state.valid){
            p["rmsd"] = state.rmsd;
        }
        perAtom.push_back(std::move(p));
    }
    return perAtom;
}

}
