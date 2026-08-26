# An Effective Memetic Algorithm for the Close-Enough Traveling Salesman Problem

This repository is the implementation in C++ for the paper [An Effective Memetic Algorithm for the Close-Enough Traveling Salesman Problem](https://www.sciencedirect.com/science/article/pii/S1568494624000401) by Zhenyu Lei, Jin-Kao Hao.

We propose an effective memetic algorithm to solve the Close-Enough Traveling Salesman Problem (CETSP) which is a variant of the well-known Traveling Salesman Problem (TSP). Experimental results on the well-known benchmark instances show that the algorithm is highly competitive with the state-of-the-art methods. We also demonstrate the usefulness of the algorithm on a real laser welding robot path planning problem.

## Overview

This vendored source tree contains:

- `cmake/` contains Gurobi discovery support;
- `include/` and `src/` contain the algorithm;
- `INTEGRATION.md` documents the custom-instance and seeded-population formats.

Benchmark data is available in the parent CETSP project's `data/` directory.

## Requirements

### Instances

We conducted experiments on the well-known benchmark instances proposed by Dr. Mennell. The instances are available at [link](https://drum.lib.umd.edu/handle/1903/9822).

And we also tested our algorithm on real-world instances provided by Prof. Béla Vizvári, one of authors of [Nedjatia et al.](https://www.jms.procedia.org/archive/CRPASE_169/CRPASE_procedia_2020_6_4_11.PDF). We provide the instances in the folder `datasets/`.

### Solver

Install [LKH](http://akira.ruc.dk/~keld/research/LKH/) to optimize visiting
sequences and [Gurobi](https://www.gurobi.com/) to optimize positions for a
fixed sequence. The CMake project recognizes Gurobi 9.5 through 13.0 library
names and has been tested with Gurobi 12.0.2.

## Run the code

Configure Gurobi through `GUROBI_ROOT` and build with a compiler compatible
with the installed Gurobi C++ library. For example, on Windows:

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 `
  -DGUROBI_ROOT=C:/gurobi1202/win64
cmake --build build --config Release
```

The original indexed-instance mode still uses defaults from `include/Defs.hpp`.
Runtime paths can instead be supplied explicitly:

```powershell
./build/Release/MA-CETSP.exe `
  --instance_file=C:/path/instance.txt `
  --seed_file=C:/path/population.seeds `
  --result_file=C:/path/result.txt `
  --lkh_executable=C:/path/LKH.exe `
  --lkh_temp_root=C:/path/lkh-temp
```

Search parameters remain available as command-line options:

```bash
./MA-CETSP -i <instance> -s <seed> -r <max-iteration> -p <population-size> -b <fitness-beta> -d <distance-threshold> -n <neighbor-size>
```

## Citation

If it is helpful for your research, please cite our paper:

```bibtex
@article{lei2024effective,
  title={An effective memetic algorithm for the close-enough traveling salesman problem},
  author={Lei, Zhenyu and Hao, Jin-Kao},
  journal={Applied Soft Computing},
  pages={111266},
  year={2024},
  publisher={Elsevier}
}
```
