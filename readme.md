# CETSP

This is an attempt at writing a fast and good algorithm for the close-enough traveling salesman problem. 

## How it works
 
Basically, it uses a heuristic based on the [pair-center algorithm](https://www.sciencedirect.com/science/article/pii/S1877750324002175). Since this is the CETSP, we also pop and re-insert points sometimes to get better results.

## Results

Currently this algorithm achieves ~10-20% worse tours than the optimal algorithms on the Mennell dataset. 

## Problems

1. It is much slower than it should be, probably due to all the adding and deleting to R-trees that I do. Theoretically it should be $O(n \times polylog(n))$, but constants don't care about your theoretical analysis. 

2. It crashes sometimes. I am trying to fix this. 