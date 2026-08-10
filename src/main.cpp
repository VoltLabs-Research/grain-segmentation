#include <volt/cli/common.h>
#include <volt/grain_segmentation_service.h>
#include <volt/plugin/option_reader.h>
#include <oneapi/tbb/global_control.h>
#include <tbb/info.h>

#include <algorithm>
#include <fstream>
#include <set>

using namespace Volt;
using namespace Volt::CLI;
using namespace Volt::Plugin;

static PluginDescriptor buildDescriptor() {
    return {
        "grain-segmentation",
        "Grain Segmentation",
        {
            {"--rmsd", "float", "RMSD threshold for PTM.", "0.1", {}, ""},
            {"--minGrainAtomCount", "int", "Minimum atoms per grain.", "100", {}, ""},
            {"--adoptOrphanAtoms", "bool", "Adopt orphan atoms.", "true", {}, ""},
            {"--handleCoherentInterfaces", "bool", "Handle coherent interfaces.", "true", {}, ""},
            {"--mergeAlgorithm", "enum", "Grain merge algorithm.", "GraphClusteringAutomatic",
             {"GraphClusteringAutomatic", "GraphClusteringManual", "MinimumSpanningTree"}, ""},
            {"--mergingThreshold", "float", "Merge threshold (used by Manual/MST modes).", "0", {}, ""},
            {"--outputBonds", "bool", "Output neighbor bonds.", "false", {}, ""},
        }
    };
}

int main(int argc, char* argv[]) {
    const PluginDescriptor descriptor = buildDescriptor();

    if (argc < 2) {
        showPluginUsage(argv[0], descriptor);
        return 1;
    }

    std::string filename, outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);

    if (auto exitCode = handleIntrospection(argv[0], descriptor, opts, filename)) {
        return *exitCode;
    }

    if (!hasOption(opts, "--threads")) {
        const int maxAvailableThreads = static_cast<int>(oneapi::tbb::info::default_concurrency());
        int physicalCores = 0;
        std::ifstream cpuinfo("/proc/cpuinfo");
        if (cpuinfo.is_open()) {
            std::set<std::pair<int, int>> physicalCoreIds;
            int fallbackCpuCores = 0;
            int physicalId = -1;
            int coreId = -1;
            std::string line;
            while (std::getline(cpuinfo, line)) {
                if (line.empty()) {
                    if (physicalId >= 0 && coreId >= 0) {
                        physicalCoreIds.emplace(physicalId, coreId);
                    }
                    physicalId = -1;
                    coreId = -1;
                    continue;
                }
                if (line.rfind("physical id", 0) == 0) {
                    physicalId = std::stoi(line.substr(line.find(':') + 1));
                } else if (line.rfind("core id", 0) == 0) {
                    coreId = std::stoi(line.substr(line.find(':') + 1));
                } else if (line.rfind("cpu cores", 0) == 0) {
                    fallbackCpuCores = std::max(fallbackCpuCores, std::stoi(line.substr(line.find(':') + 1)));
                }
            }
            if (physicalId >= 0 && coreId >= 0) {
                physicalCoreIds.emplace(physicalId, coreId);
            }
            physicalCores = !physicalCoreIds.empty()
                ? static_cast<int>(physicalCoreIds.size())
                : fallbackCpuCores;
        }
        int defaultThreads = maxAvailableThreads;
        if (physicalCores > 0) {
            defaultThreads = std::min(maxAvailableThreads, physicalCores);
        }
        opts["--threads"] = std::to_string(std::max(1, defaultThreads));
    }

    const int requestedThreads = getInt(opts, "--threads");
    oneapi::tbb::global_control parallelControl(
        oneapi::tbb::global_control::max_allowed_parallelism,
        static_cast<std::size_t>(std::max(1, requestedThreads))
    );
    initLogging("grain-segmentation");
    spdlog::info("Using {} threads (OneTBB)", requestedThreads);
    
    LammpsParser::Frame frame;
    if (!parseFrame(filename, frame)) return 1;
    
    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);
    
    const OptionReader options(descriptor, opts);

    bool adoptOrphanAtoms = options.boolean("--adoptOrphanAtoms");
    int minGrainAtomCount = options.integer("--minGrainAtomCount");
    bool handleCoherentInterfaces = options.boolean("--handleCoherentInterfaces");
    bool outputBonds = options.boolean("--outputBonds");

    const std::string mergeAlgorithmStr = options.text("--mergeAlgorithm");
    MergeAlgorithm mergeAlgorithm = MergeAlgorithm::GraphClusteringAutomatic;
    if(mergeAlgorithmStr == "GraphClusteringManual"){
        mergeAlgorithm = MergeAlgorithm::GraphClusteringManual;
    }else if(mergeAlgorithmStr == "MinimumSpanningTree"){
        mergeAlgorithm = MergeAlgorithm::MinimumSpanningTree;
    }else if(mergeAlgorithmStr != "GraphClusteringAutomatic"){
        spdlog::warn("Unknown mergeAlgorithm '{}', defaulting to GraphClusteringAutomatic", mergeAlgorithmStr);
    }
    double mergingThreshold = options.number("--mergingThreshold");

    spdlog::info("Grain segmentation parameters:");
    spdlog::info("  - adoptOrphanAtoms: {}", adoptOrphanAtoms);
    spdlog::info("  - minGrainAtomCount: {}", minGrainAtomCount);
    spdlog::info("  - handleCoherentInterfaces: {}", handleCoherentInterfaces);
    spdlog::info("  - mergeAlgorithm: {}", mergeAlgorithmStr);
    spdlog::info("  - mergingThreshold: {}", mergingThreshold);
    spdlog::info("  - outputBonds: {}", outputBonds);

    GrainSegmentationService analyzer;
    analyzer.setRMSD(options.number("--rmsd"));
    analyzer.setParameters(
        adoptOrphanAtoms,
        minGrainAtomCount,
        handleCoherentInterfaces,
        outputBonds,
        mergeAlgorithm,
        mergingThreshold
    );
    
    spdlog::info("Starting grain segmentation...");
    json result = analyzer.compute(frame, outputBase);
    
    if (result.value("is_failed", false)) {
        spdlog::error("Analysis failed: {}", result.value("error", "Unknown error"));
        return 1;
    }
    return 0;
}
