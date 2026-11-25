#include <volt/plugin/plugin_entry.h>
#include <volt/grain_segmentation_service.h>

using namespace Volt;
using namespace Volt::Plugin;
using S = GrainSegmentationService;

static const std::vector<OptionBinding<S>> bindings = {
    opt("--rmsd", "RMSD threshold for PTM.", 0.1f, &S::setRMSD),
    opt("--minGrainAtomCount", "Minimum atoms per grain.", 100, &S::setMinGrainAtomCount),
    opt("--adoptOrphanAtoms", "Adopt orphan atoms.", true, &S::setAdoptOrphanAtoms),
    opt("--handleCoherentInterfaces", "Handle coherent interfaces.", true, &S::setHandleCoherentInterfaces),
    opt({"--mergeAlgorithm", "enum", "Grain merge algorithm.", "GraphClusteringAutomatic",
         {"GraphClusteringAutomatic", "GraphClusteringManual", "MinimumSpanningTree"}, ""},
        std::function<void(S&, const OptsMap&)>([](S& s, const OptsMap& opts) {
            const std::string name = CLI::getString(opts, "--mergeAlgorithm", "GraphClusteringAutomatic");
            MergeAlgorithm algorithm = MergeAlgorithm::GraphClusteringAutomatic;
            if (name == "GraphClusteringManual") {
                algorithm = MergeAlgorithm::GraphClusteringManual;
            } else if (name == "MinimumSpanningTree") {
                algorithm = MergeAlgorithm::MinimumSpanningTree;
            } else if (name != "GraphClusteringAutomatic") {
                spdlog::warn("Unknown mergeAlgorithm '{}', defaulting to GraphClusteringAutomatic", name);
            }
            s.setMergeAlgorithm(algorithm);
        })),
    opt({"--mergingThreshold", "float", "Merge threshold (used by Manual/MST modes).", "0", {}, ""},
        std::function<void(S&, const OptsMap&)>([](S& s, const OptsMap& opts) {
            s.setMergingThreshold(CLI::getDouble(opts, "--mergingThreshold", 0.0));
        })),
    opt("--outputBonds", "Output neighbor bonds.", false, &S::setOutputBonds),
};

VOLT_PLUGIN_MAIN((PluginDescriptor{
        .name = "grain-segmentation",
        .description = "Grain Segmentation",
        .options = optionsMeta(bindings),
        .defaultThreads = DefaultThreads::PhysicalCores,
    }),
    serviceRun(bindings))
