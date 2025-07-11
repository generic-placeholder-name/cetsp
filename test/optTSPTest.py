import math
import numpy as np
import matplotlib.pyplot as plt
from functools import lru_cache

# 1) Read 2D points from file (or generate random points if file not found)
file_path = '../data/cmu_output/rotatingDiamonds5.txt'
points = []
try:
    with open(file_path) as f:
        for line in f:
            x, y = map(float, line.strip().split())
            points.append((x, y))
    print(f"Loaded {len(points)} points from '{file_path}'.")
except FileNotFoundError:
    n = 50
    points = [(np.random.rand(), np.random.rand()) for _ in range(n)]
    print(f"File not found. Generated {n} random points for demonstration.")

n = len(points)
tour = list(range(n))

# 2) Distance function
def dist(i, j):
    return math.hypot(points[i][0] - points[j][0],
                      points[i][1] - points[j][1])

# 3) Compute original tour length (including last→first)
old_distance = sum(dist(i, (i + 1) % n) for i in range(n))

# 4) Compute k = round(log2(n) - log2(log2(n))), bounded [1, n-1]
k = int(round(math.log2(n) - math.log2(math.log2(n))))
k = max(1, min(k, n-1))
print(f"Using k = {k}.")

# 5) Identify the k largest edges in the original tour
edges = [(i, (i + 1) % n, dist(i, (i + 1) % n)) for i in range(n)]
edges_sorted = sorted(edges, key=lambda x: x[2], reverse=True)
to_remove = set((i, j) for i, j, _ in edges_sorted[:k])

# 6) Break tour into segments by removing those k edges
visited = [False] * n
segments = []
for i in range(n):
    if not visited[i]:
        seg = [i]
        visited[i] = True
        cur = i
        while True:
            nxt = (cur + 1) % n
            if (cur, nxt) in to_remove or visited[nxt]:
                break
            seg.append(nxt)
            visited[nxt] = True
            cur = nxt
        segments.append(seg)

m = len(segments)
print(f"Created {m} segments after removing edges.")

# 7) Precompute segment endpoints (start, end) for both orientations
endpoints = [(points[seg[0]], points[seg[-1]]) for seg in segments]

# 8) Precompute cost between every pair of segments for both orientations
cost = {}
for i in range(m):
    cost[i] = {}
    for j in range(m):
        if i == j: continue
        cost[i][j] = {}
        for oi in (0, 1):
            for oj in (0, 1):
                # if oi==0, i goes start→end; if oi==1, it's reversed end→start
                end_i = endpoints[i][1] if oi == 0 else endpoints[i][0]
                # if oj==0, j goes start→end; if oj==1, reversed
                start_j = endpoints[j][0] if oj == 0 else endpoints[j][1]
                cost[i][j][(oi, oj)] = math.hypot(
                    end_i[0] - start_j[0],
                    end_i[1] - start_j[1]
                )

# 9) Precompute closing‐edge cost from every segment back to segment 0 (orientation fixed to 0)
close_cost = {}
for last in range(m):
    for ori in (0, 1):
        end_last = endpoints[last][1] if ori == 0 else endpoints[last][0]
        start_first = endpoints[0][0]  # segment 0, orientation 0
        close_cost[(last, ori)] = math.hypot(
            end_last[0] - start_first[0],
            end_last[1] - start_first[1]
        )

# 10) DP over subsets of segments, carrying (mask, last_seg, last_ori)
@lru_cache(None)
def dp(mask, last, ori):
    # when all segments are used, pay the closing‐edge cost
    if mask == (1 << m) - 1:
        return close_cost[(last, ori)]
    best = float('inf')
    for j in range(m):
        if mask & (1 << j):
            continue
        for oj in (0, 1):
            c = cost[last][j][(ori, oj)] + dp(mask | (1 << j), j, oj)
            if c < best:
                best = c
    return best

# 11) Run DP starting from segment 0 with orientation 0
total_new_cost = dp(1 << 0, 0, 0)

# 12) Reconstruct optimal segment order & orientations
mask, last, ori = (1 << 0), 0, 0
order = [(0, 0)]
while mask != (1 << m) - 1:
    for j in range(m):
        if mask & (1 << j):
            continue
        for oj in (0, 1):
            if abs(cost[last][j][(ori, oj)]
                   + dp(mask | (1 << j), j, oj)
                   - dp(mask, last, ori)) < 1e-6:
                order.append((j, oj))
                mask |= (1 << j)
                last, ori = j, oj
                break
        else:
            continue
        break

# 13) Build the new tour index list (and close it)
new_tour = []
for seg_idx, orientation in order:
    seg = segments[seg_idx][:]
    if orientation == 1:
        seg.reverse()
    new_tour.extend(seg)
new_tour.append(new_tour[0])
new_distance = sum(dist(new_tour[i], new_tour[i + 1]) for i in range(len(new_tour) - 1))

old_x = [points[i][0] for i in tour] + [points[tour[0]][0]]
old_y = [points[i][1] for i in tour] + [points[tour[0]][1]]
plt.plot(old_x, old_y)
plt.title('Old Tour')

plt.figure()
new_x = [points[i][0] for i in new_tour]
new_y = [points[i][1] for i in new_tour]
plt.plot(new_x, new_y)
plt.title('New Tour')

plt.figure()
plt.bar(['Old', 'New'], [old_distance, new_distance])
plt.title('Tour Length Comparison')

plt.show()

# 15) Print results
print(f"Old tour length: {old_distance:.4f}")
print(f"New tour length: {new_distance:.4f}")
