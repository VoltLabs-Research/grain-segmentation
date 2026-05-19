# GrainSegmentation

`GrainSegmentation` segments grains from the current frame using PTM-backed local structure information.

## One-Command Install

```bash
curl -sSL https://raw.githubusercontent.com/VoltLabs-Research/CoreToolkit/main/scripts/install-plugin.sh | bash -s -- GrainSegmentation
```

## Build from source

Requires [Conan 2.x](https://docs.conan.io/2/installation.html), CMake 3.20+, and a C++23 compiler (GCC 14+ or Clang 17+).

### Prerequisites

The following Conan packages must be available in your local cache:

- `coretoolkit/1.0.0` (from the `CoreToolkit` repository)
- `structure-identification/1.0.0` (from the `StructureIdentification` repository)
- `polyhedral-template-matching/1.0.0` (from the `PolyhedralTemplateMatching` repository)

For each dependency, clone its repository and create the package:

```bash
conan create <path-to-dependency-repo> --build=missing -o "hwloc/*:shared=True"
```

### Build

From the root of this repository:

```bash
conan install . -of build --build=missing -o "hwloc/*:shared=True"
cmake --preset conan-release
cmake --build build/build/Release -j
```

### Run

```bash
./build/build/Release/grain-segmentation --help
```

### Package as Conan recipe

To make this plugin available as a Conan package for other projects:

```bash
conan create . --build=missing -o "hwloc/*:shared=True"
```

## CLI

Usage:

```bash
grain-segmentation <lammps_file> [output_base] [options]
```

### Arguments

| Argument | Required | Description | Default |
| --- | --- | --- | --- |
| `<lammps_file>` | Yes | Input LAMMPS dump file. | |
| `[output_base]` | No | Base path for output files. | derived from input |
| `--rmsd <float>` | No | RMSD threshold for PTM. | `0.1` |
| `--minGrainAtomCount <int>` | No | Minimum atoms per grain. | `100` |
| `--adoptOrphanAtoms <true\|false>` | No | Adopt orphan atoms into neighboring grains. | `true` |
| `--handleCoherentInterfaces <true\|false>` | No | Handle coherent interfaces specially. | `true` |
| `--outputBonds` | No | Export neighbor bonds. | `false` |
| `--threads <int>` | No | Maximum worker threads. | auto |
| `--help` | No | Print CLI help. | |
