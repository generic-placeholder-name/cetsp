# CETSP

This is an attempt at writing a fast and good algorithm for the close-enough traveling salesman problem. 

## How it works

Basically, it uses a heuristic based on the [pair-center algorithm](https://www.sciencedirect.com/science/article/pii/S1877750324002175). We iteratively 'merge' two circles into a 'representative' circle until there is one left. Then, we 'pop' one circle into two, maintaining a tour, until we get the original circles again. 

Furthermore, every once in a while, a node will be deleted, and every circle corresponding to that node will be re-inserted into the tree. Also, occasionally, a node will have its position optimized based on its two neighboring points.

## Requirements

You must have Boost installed. Abseil is also used for fast hash maps; you can use the STL ones instead, but it is slower. 

## Building (Windows, MinGW; UNIX/POSIX similar)

1. Open a terminal (PowerShell, Command Prompt, or UNIX shell).
2. Make a build directory (if it doesn’t exist) and navigate into it:

```bash
mkdir build   # Windows: mkdir build
cd build
```

3. Configure the project with CMake using **MinGW**:

```bash
cmake .. -G "MinGW Makefiles"
```

4. Build the project:

```bash
cmake --build .
```

5. Go back to the base directory:

```bash
cd ..
```

**Notes:**

* On UNIX/POSIX, you can replace `MinGW Makefiles` with `Unix Makefiles` or use `Ninja` (`-G "Ninja"`).
* You can also pass options to CMake:

```bash
cmake .. -DUSE_ABSEIL_HASH_SET=ON -DDEBUG=OFF -G "MinGW Makefiles"
```

* `USE_ABSEIL_HASH_SET=ON` → uses Abseil flat_hash_set (`ON` by default).
* `DEBUG=ON` → enables debug logging (`DBG()` macros) (`OFF` by default). (Warning: the amount of debug data is massive. This is mostly used to help debug code when it's not working.)

## Running

* Place your `settings.txt` file in the **base directory** (same level as the `build` folder).
* Run the executable:

```bash
./build/CETSP_project.exe
```

* On Windows, you can also double-click the `.exe`.

## Results

Currently, this algorithm achieves approximately 0–5% worse tours than SOTA on the Mennell dataset. Visualized tours can be found in [`data/output_images`](./data/output_images). 

## Problems

The TSP order itself is sometimes suboptimal. Fixing some of this probably requires getting a better TSP algorithm, but those incur quadratic or larger runtime. Others probably demand a better placement or popping heuristic.