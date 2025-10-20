import os
import random
import math

def generate_structured_cetsp(output_dir="timing_structured_input",
                              exps=range(7, 21),
                              jitter=0.1,
                              radius_min=0.2,
                              radius_max=0.5,
                              rng_seed=None):
    """
    Generate structured CETSP inputs:
      - For each n = 2^exp, place n centers roughly on a sqrt(n) x sqrt(n) grid,
        with small jitter in [-jitter, +jitter].
      - Radii are sampled uniformly in [radius_min, radius_max].
    """
    if rng_seed is not None:
        random.seed(rng_seed)

    os.makedirs(output_dir, exist_ok=True)

    for exp in exps:
        n = 2 ** exp
        side = int(math.ceil(math.sqrt(n))) - 1 # so that we have some extra random points to play with
        filename = os.path.join(output_dir, f"test_n{n}.txt")
        centers = []

        # Fill grid row-major until we have n points
        count = 0
        for i in range(side):
            for j in range(side):
                if count >= n:
                    break
                # base grid position (i, j) scaled so average spacing ~1
                # If you want spacing > 1, multiply i,j by spacing
                x = i + random.uniform(-jitter, jitter)
                y = j + random.uniform(-jitter, jitter)
                r = random.uniform(radius_min, radius_max)
                centers.append((x, y, r))
                count += 1
            if count >= n:
                break

        # Sanity: if side^2 > n we only used first n cells; otherwise we've filled exactly n
        with open(filename, "w") as f:
            for (x, y, r) in centers:
                f.write(f"{x:.6f} {y:.6f} {r:.6f}\n")

        print(f"Generated {filename} with {n} circles (grid side={side}).")

if __name__ == "__main__":
    generate_structured_cetsp(rng_seed=12345)
