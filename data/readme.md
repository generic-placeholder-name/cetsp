# Dataset

This repository contains the **CMU/Mennell dataset** for the CETSP.

- The dataset located at `./cmu` is the publicly available version from the [UMD repository](https://drum.lib.umd.edu/items/a532bfe4-9fe2-40ab-b562-f5851ddd87d9).
- Note: This dataset is **not identical** to the one used to benchmark previous papers. I suspect that the dataset was changed some time ago. 

To process the dataset into a comparable format, use the provided Python scripts in this directory.

The `benchmark` directory, from Lei and Hao (2024), represents SOTA solutions.


## Results

### Table 1: Comparison of benchmark vs output values for all instances

| Instance | N | OR% | Benchmark Value | Output Value <sup>[1]</sup> | Gap% |
| --- | --- | --- | --- | --- | --- |
| bonus1000 <sup>[2]</sup> | 1000 | 12.26 | 384.365 | 378.622 | -1.49 |
| bonus1000rdmRad | 1000 | 6.13 | 917.610 | 937.247 | 2.14 |
| bubbles1 | 36 | 11.11 | 349.135 | 349.255 | 0.03 |
| bubbles2 | 76 | 9.09 | 428.279 | 428.368 | 0.02 |
| bubbles3 | 126 | 7.69 | 529.955 | 531.243 | 0.24 |
| bubbles4 | 184 | 6.67 | 802.974 | 809.216 | 0.78 |
| bubbles5 | 250 | 5.88 | 1035.320 | 1055.488 | 1.95 |
| bubbles6 | 324 | 5.26 | 1220.070 | 1301.459 | 6.67 |
| bubbles7 | 406 | 4.76 | 1575.040 | 1647.837 | 4.62 |
| bubbles8 | 496 | 4.35 | 1881.930 | 1995.234 | 6.02 |
| bubbles9 | 594 | 4.00 | 2148.400 | 2331.431 | 8.52 |
| chaoSingleDep | 200 | 2.84 | 1039.610 | 1042.864 | 0.31 |
| concentricCircles1 | 16 | 15.00 | 53.158 | 53.201 | 0.08 |
| concentricCircles2 | 36 | 7.50 | 153.132 | 153.870 | 0.48 |
| concentricCircles3 | 60 | 5.00 | 270.007 | 271.744 | 0.64 |
| concentricCircles4 | 104 | 3.75 | 451.870 | 460.079 | 1.82 |
| concentricCircles5 | 148 | 3.00 | 632.977 | 639.080 | 0.96 |
| d493_or10 | 493 | 10.00 | 100.721 | 100.959 | 0.24 |
| d493_or2 | 493 | 2.00 | 199.015 | 200.798 | 0.90 |
| d493_or30 | 493 | 30.00 | 69.758 | 69.758 | 0.00 |
| d493rdmRad | 493 | 13.44 | 134.226 | 134.621 | 0.29 |
| dsj1000_or10 | 1000 | 10.00 | 373.759 | 376.445 | 0.72 |
| dsj1000_or2 | 1000 | 2.00 | 909.502 | 922.440 | 1.42 |
| dsj1000_or30 | 1000 | 30.00 | 199.948 | 200.037 | 0.04 |
| dsj1000rdmRad | 1000 | 12.47 | 624.604 | 626.148 | 0.25 |
| kroD100_or10 | 100 | 10.00 | 89.668 | 89.781 | 0.13 |
| kroD100_or2 | 100 | 2.00 | 159.037 | 159.603 | 0.36 |
| kroD100_or30 | 100 | 30.00 | 58.541 | 58.576 | 0.06 |
| kroD100rdmRad | 100 | 4.03 | 141.829 | 142.199 | 0.26 |
| lin318_or10 | 318 | 10.00 | 1394.630 | 1404.343 | 0.70 |
| lin318_or2 | 318 | 2.00 | 2816.580 | 2833.410 | 0.60 |
| lin318_or30 | 318 | 30.00 | 765.964 | 765.998 | 0.00 |
| lin318rdmRad | 318 | 10.05 | 2047.110 | 2052.486 | 0.26 |
| pcb442_or10 | 442 | 10.00 | 142.573 | 143.302 | 0.51 |
| pcb442_or2 | 442 | 2.00 | 319.846 | 332.701 | 4.02 |
| pcb442_or30 | 442 | 30.00 | 83.537 | 83.538 | 0.00 |
| pcb442rdmRad | 442 | 9.39 | 219.219 | 220.280 | 0.48 |
| rat195_or10 | 195 | 10.00 | 67.991 | 68.080 | 0.13 |
| rat195_or2 | 195 | 2.00 | 157.970 | 160.495 | 1.60 |
| rat195_or30 | 195 | 30.00 | 45.702 | 45.702 | 0.00 |
| rat195rdmRad | 195 | 42.60 | 68.224 | 68.252 | 0.04 |
| rd400_or10 | 400 | 10.00 | 452.826 | 453.943 | 0.25 |
| rd400_or2 | 400 | 2.00 | 1018.460 | 1040.517 | 2.17 |
| rd400_or30 | 400 | 30.00 | 224.839 | 224.910 | 0.03 |
| rd400rdmRad | 400 | 0.99 | 1238.280 | 1279.515 | 3.33 |
| rotatingDiamonds1 | 20 | 20.00 | 32.389 | 32.431 | 0.13 |
| rotatingDiamonds2 | 60 | 5.00 | 140.477 | 140.636 | 0.11 |
| rotatingDiamonds3 | 180 | 3.33 | 380.882 | 381.475 | 0.16 |
| rotatingDiamonds4 | 320 | 1.43 | 770.661 | 790.601 | 2.59 |
| rotatingDiamonds5 | 680 | 1.11 | 1510.750 | 1554.064 | 2.87 |
| team1_100 | 100 | 9.33 | 307.337 | 308.339 | 0.33 |
| team1_100rdmRad | 100 | 7.69 | 388.537 | 399.272 | 2.76 |
| team2_200 | 200 | 20.06 | 246.683 | 247.401 | 0.29 |
| team2_200rdmRad | 200 | 5.19 | 613.659 | 618.001 | 0.71 |
| team3_300 | 300 | 7.02 | 461.889 | 462.960 | 0.23 |
| team3_300rdmRad | 300 | 23.70 | 378.087 | 379.137 | 0.28 |
| team4_400 | 400 | 5.01 | 669.910 | 674.344 | 0.66 |
| team4_400rdmRad | 400 | 2.05 | 984.240 | 1013.142 | 2.94 |
| team5_499 | 499 | 2.00 | 693.798 | 704.900 | 1.60 |
| team5_499rdmRad | 499 | 20.14 | 446.191 | 446.525 | 0.07 |
| team6_500 | 500 | 27.06 | 225.216 | 225.451 | 0.10 |
| team6_500rdmRad | 500 | 10.02 | 620.886 | 623.678 | 0.45 |

**Notes:**

* [1] Output values obtained using our CETSP solver over 1000 runs; the best value is reported.
* [2] The `bonus1000` test case was ran 10000 times to improve the best value.

---

### Table 2: Runtime comparison: benchmark vs our run (adjusted)

| Instance | Benchmark Time (s) | Our time (s) | Adjusted time (s) <sup>[1]</sup> | Adj / Bench <sup>[2]</sup> |
| --- | --- | --- | --- | --- |
| bonus1000 | 2016.81 | 0.300 | 3.127 | 0.0016 |
| bonus1000rdmRad | 890.66 | 0.239 | 2.491 | 0.0028 |
| bubbles1 | 31.12 | 0.014 | 0.146 | 0.0047 |
| bubbles2 | 41.73 | 0.025 | 0.261 | 0.0062 |
| bubbles3 | 193.23 | 0.051 | 0.532 | 0.0028 |
| bubbles4 | 173.63 | 0.077 | 0.803 | 0.0046 |
| bubbles5 | 247.52 | 0.125 | 1.303 | 0.0053 |
| bubbles6 | 358.74 | 0.192 | 2.001 | 0.0056 |
| bubbles7 | 650.29 | 0.228 | 2.377 | 0.0037 |
| bubbles8 | 790.10 | 0.292 | 3.044 | 0.0039 |
| bubbles9 | 893.73 | 0.443 | 4.618 | 0.0052 |
| chaoSingleDep | 76.44 | 0.059 | 0.615 | 0.0080 |
| concentricCircles1 | 29.31 | 0.002 | 0.021 | 0.0007 |
| concentricCircles2 | 139.32 | 0.011 | 0.115 | 0.0008 |
| concentricCircles3 | 67.01 | 0.023 | 0.240 | 0.0036 |
| concentricCircles4 | 193.29 | 0.049 | 0.511 | 0.0026 |
| concentricCircles5 | 208.92 | 0.082 | 0.855 | 0.0041 |
| d493_or10 | 352.52 | 0.033 | 0.344 | 0.0010 |
| d493_or2 | 689.44 | 0.140 | 1.459 | 0.0021 |
| d493_or30 | 116.74 | 0.267 | 2.783 | 0.0238 |
| d493rdmRad | 71.85 | 0.189 | 1.970 | 0.0274 |
| dsj1000_or10 | 942.49 | 0.089 | 0.928 | 0.0010 |
| dsj1000_or2 | 2578.76 | 0.240 | 2.502 | 0.0010 |
| dsj1000_or30 | 249.61 | 0.597 | 6.223 | 0.0249 |
| dsj1000rdmRad | 176.72 | 0.239 | 2.491 | 0.0141 |
| kroD100_or10 | 50.76 | 0.030 | 0.313 | 0.0062 |
| kroD100_or2 | 158.32 | 0.020 | 0.208 | 0.0013 |
| kroD100_or30 | 39.25 | 0.038 | 0.396 | 0.0101 |
| kroD100rdmRad | 158.19 | 0.009 | 0.094 | 0.0006 |
| lin318_or10 | 193.10 | 0.047 | 0.490 | 0.0025 |
| lin318_or2 | 327.54 | 0.094 | 0.981 | 0.0030 |
| lin318_or30 | 71.47 | 0.199 | 2.076 | 0.0290 |
| lin318rdmRad | 66.44 | 0.059 | 0.615 | 0.0093 |
| pcb442_or10 | 467.38 | 0.082 | 0.855 | 0.0018 |
| pcb442_or2 | 1334.95 | 0.033 | 0.344 | 0.0003 |
| pcb442_or30 | 117.86 | 0.140 | 1.459 | 0.0124 |
| pcb442rdmRad | 160.56 | 0.267 | 2.783 | 0.0173 |
| rat195_or10 | 170.90 | 0.189 | 1.970 | 0.0115 |
| rat195_or2 | 322.37 | 0.089 | 0.928 | 0.0029 |
| rat195_or30 | 50.19 | 0.240 | 2.502 | 0.0499 |
| rat195rdmRad | 32.17 | 0.597 | 6.223 | 0.1934 |
| rd400_or10 | 540.53 | 0.239 | 2.491 | 0.0046 |
| rd400_or2 | 450.63 | 0.030 | 0.313 | 0.0007 |
| rd400_or30 | 91.67 | 0.020 | 0.208 | 0.0023 |
| rd400rdmRad | 635.69 | 0.038 | 0.396 | 0.0006 |
| rotatingDiamonds1 | 33.70 | 0.059 | 0.615 | 0.0182 |
| rotatingDiamonds2 | 46.81 | 0.082 | 0.855 | 0.0183 |
| rotatingDiamonds3 | 87.69 | 0.105 | 1.094 | 0.0125 |
| rotatingDiamonds4 | 229.92 | 0.302 | 3.149 | 0.0137 |
| rotatingDiamonds5 | 597.20 | 0.073 | 0.761 | 0.0013 |
| team1_100 | 54.29 | 0.009 | 0.094 | 0.0017 |
| team1_100rdmRad | 38.02 | 0.050 | 0.521 | 0.0137 |
| team2_200 | 181.62 | 0.122 | 1.270 | 0.0070 |
| team2_200rdmRad | 93.86 | 0.022 | 0.230 | 0.0024 |
| team3_300 | 435.69 | 0.287 | 2.992 | 0.0069 |
| team3_300rdmRad | 48.46 | 0.133 | 1.389 | 0.0287 |
| team4_400 | 648.08 | 0.273 | 2.844 | 0.0044 |
| team4_400rdmRad | 734.98 | 0.055 | 0.574 | 0.0008 |
| team5_499 | 561.82 | 0.003 | 0.031 | 0.0001 |
| team5_499rdmRad | 51.36 | 0.026 | 0.271 | 0.0053 |
| team6_500 | 208.22 | 0.114 | 1.188 | 0.0057 |
| team6_500rdmRad | 261.21 | 0.229 | 2.384 | 0.0091 |

**Notes:**

* [1] Our runtimes are measured on an AMD Ryzen 9 7940HS and reported in seconds. Times are scaled to the benchmark CPU (AMD Opteron 4184) using a [PassMark](https://passmark.com)-based factor of 10.424 (29938 / 2872), following Lei & Hao (2024).
* [2] Adj / Bench = Adjusted time / Benchmark time.