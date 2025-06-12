# CETSP

This is an attempt at writing a fast and good algorithm for the close-enough traveling salesman problem. 

## How it works
 
Basically, it uses a heuristic based on the [pair-center algorithm](https://www.sciencedirect.com/science/article/pii/S1877750324002175). We iteratively 'merge' two circles into a 'representative' circle until there is one left. Then, we 'pop' one circle into two, maintaining a tour, until we get the original circles again. 

Furthermore, every once in a while, a node will be deleted, and every circle corresponding to that node will be re-inserted into the tree. Also, occasionally, a node will have its position optimized based on its two neighboring points.

## Results

Currently, this algorithm achieves approximately 0–5% worse tours than SOTA on the Mennell dataset. Visualized tours can be found in [`data/output_images`](./data/output_images). 

## Problems

1. This is supposed to be $O(n \times polylog(n))$, but since I use `std::vector` in my `TourNodes`, it can be $O(n^2)$. Fixing this is simple and will be done soonish.
2. The TSP order itself is sometimes suboptimal. Fixing this probably requires getting a better TSP algorithm, but those incur quadratic or larger runtime.