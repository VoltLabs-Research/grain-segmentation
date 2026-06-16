# Grain Segmentation

Segments crystalline grains from local PTM structure, merging regions with similar orientations.

## Install

```bash
vpm install @voltlabs/grain-segmentation
```

## CLI

```bash
grain-segmentation <input_dump> [output_base] [options]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `<input_dump>` | yes | — | Input LAMMPS dump. |
| `[output_base]` | no | derived from input | Base path for output files. |
| `--rmsd <float>` | no | `0.1` | RMSD threshold for PTM. |
| `--minGrainAtomCount <int>` | no | `100` | Minimum atoms per grain. |
| `--adoptOrphanAtoms <true\|false>` | no | `true` | Adopt orphan atoms into neighboring grains. |
| `--handleCoherentInterfaces <true\|false>` | no | `true` | Handle coherent interfaces specially. |
| `--mergeAlgorithm <name>` | no | `GraphClusteringAutomatic` | Merge algorithm: `GraphClusteringAutomatic`, `GraphClusteringManual`, or `MinimumSpanningTree`. |
| `--mergingThreshold <float>` | no | `0` | Merge threshold (used by Manual/MST modes). |
| `--outputBonds` | no | `false` | Export neighbor bonds. |
| `--threads <int>` | no | auto | Max worker threads. |

## Exports

| Output file | Exposure | Exporter → artifact |
|---|---|---|
| `{output_base}_grains.parquet` | Grain Segmentation | — (listing-only) |
| `{output_base}_atoms.parquet` | Grain Model | AtomisticExporter → glb |

---

Full input contract and examples: https://docs.voltcloud.dev/docs/plugins/grain-segmentation
