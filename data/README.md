# Mennell benchmark data

This directory contains the CETSP instances used in the repository's Mennell
benchmark comparisons.

```text
mennell/original/         public `.cetsp` files from the UMD archive
mennell/instances/        normalized `x y radius` input files
mennell/reference_tours/  tours and reported values distributed by Lei and Hao
```

The public UMD data are not necessarily byte-for-byte identical to the private
inputs used in every earlier publication. Comparisons in this repository use
the normalized instances checked in here. Dataset conversion and visualization
utilities are in [`tools/data`](../tools/data); measurements and generated tours
are kept under [`results`](../results), not alongside the inputs.

The source archive is available from the
[University of Maryland repository](https://drum.lib.umd.edu/items/a532bfe4-9fe2-40ab-b562-f5851ddd87d9).
