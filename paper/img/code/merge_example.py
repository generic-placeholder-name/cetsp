import os
import math
import random
import zipfile
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# ------------------------
# Geometry / proxy logic
# ------------------------
def make_combined_circle(p1, r1, p2, r2, deterministic=True):
    """Return a proxy circle for two input circles (same geometry as original).
       If deterministic True, choose a deterministic radius inside allowed range."""
    x1, y1 = p1
    x2, y2 = p2
    dx, dy = x2 - x1, y2 - y1
    d = math.hypot(dx, dy)
    # one inside the other
    if d + min(r1, r2) <= max(r1, r2):
        return (p1, r1) if r1 >= r2 else (p2, r2)
    ux, uy = dx / d, dy / d
    xEdge1, yEdge1 = x1 + ux * r1, y1 + uy * r1
    xEdge2, yEdge2 = x2 - ux * r2, y2 - uy * r2
    center = ((xEdge1 + xEdge2) / 2.0, (yEdge1 + yEdge2) / 2.0)
    # non-overlapping -> proxy is a point
    if d >= r1 + r2:
        return (center, 0.0)
    overlap_depth = (r1 + r2 - d) / 2.0
    a = (r1**2 - r2**2 + d**2) / (2*d)
    h = math.sqrt(max(0.0, r1**2 - a**2))
    if deterministic:
        radius = (overlap_depth + h) / 2.0
    else:
        radius = random.uniform(overlap_depth, h)
    return (center, radius)

def circle_distance(c1, c2):
    (p1, r1), (p2, r2) = c1, c2
    return math.dist(p1, p2) - r1 - r2

# ------------------------
# Plotting helpers
# ------------------------
def plot_circles_merge(ax, circles, tour_pts=None, highlight=None, proxy=None, step=None):
    ax.clear()
    ax.set_aspect('equal', 'box')
    ax.set_xlim(-7, 7)
    ax.set_ylim(-7, 7)
    ax.set_title(f"Step {step}" if step is not None else "", fontsize=12)

    # draw circles
    for i, (p, r) in enumerate(circles):
        color = 'gray'
        lw = 1.0
        if highlight and i in highlight:
            color = 'tab:red'
            lw = 2.5
        ax.add_patch(Circle(p, r, fill=False, color=color, lw=lw))
        ax.plot(p[0], p[1], 'o', color=color, markersize=4, alpha=0.9)

    # draw proxy circle (special)
    if proxy:
        (p_proxy, r_proxy) = proxy
        ax.add_patch(Circle(p_proxy, r_proxy, fill=False, color='tab:blue', lw=2.0, linestyle='--'))
        ax.plot(p_proxy[0], p_proxy[1], 'o', color='tab:blue', markersize=5)

    # draw tour (closed polyline) if provided
    if tour_pts:
        xs = [p[0] for p in tour_pts]
        ys = [p[1] for p in tour_pts]
        if len(xs) > 0:
            xs_closed = xs + [xs[0]]
            ys_closed = ys + [ys[0]]
            ax.plot(xs_closed, ys_closed, '-', linewidth=2, alpha=0.85)
            ax.plot(xs, ys, 'o', markersize=6)

    plt.tight_layout()

def plot_circles(ax,
                 circles,
                 tour_pts=None,
                 removed_circles=None,
                 removed_pts=None,
                 inserted_circles=None,
                 inserted_pts=None,
                 step=None):

    ax.clear()
    ax.set_aspect('equal', 'box')
    ax.set_xlim(-7, 7)
    ax.set_ylim(-7, 7)
    ax.set_title(f"Step {step}" if step is not None else "", fontsize=12)

    # --- existing circles ---
    for (p, r) in circles:
        ax.add_patch(Circle(p, r, fill=False, color='gray', lw=1.2))
        ax.plot(p[0], p[1], 'o', color='gray', markersize=4)

    # --- removed circles (dotted gray) ---
    if removed_circles:
        for (p, r) in removed_circles:
            ax.add_patch(Circle(p, r, fill=False, color='gray',
                                lw=1.2, linestyle=':'))
            ax.plot(p[0], p[1], 'o', color='tab:blue', markersize=4, alpha=0.6)

    # --- inserted circles (red) ---
    if inserted_circles:
        for (p, r) in inserted_circles:
            ax.add_patch(Circle(p, r, fill=False, color='tab:red', lw=2.2))
            ax.plot(p[0], p[1], 'o', color='tab:red', markersize=5)

    # --- tour ---
    if tour_pts and len(tour_pts) > 0:
        xs = [p[0] for p in tour_pts]
        ys = [p[1] for p in tour_pts]
        xs_c = xs + [xs[0]]
        ys_c = ys + [ys[0]]
        ax.plot(xs_c, ys_c, '-', lw=2, color='black')
        ax.plot(xs, ys, 'o', color='black', markersize=6)

    # --- removed tour points ---
    if removed_pts:
        for p in removed_pts:
            ax.plot(p[0], p[1], 'o', color='gray',
                    markersize=6, alpha=0.6)

    # --- inserted tour points ---
    if inserted_pts:
        for p in inserted_pts:
            ax.plot(p[0], p[1], 'o', color='tab:red',
                    markersize=7)

    plt.tight_layout()

# ------------------------
# Merge simulation (records history)
# ------------------------
def simulate_merge(n=9, seed=0, merge_save_dir="../merge_frames", deterministic_proxy=True):
    random.seed(seed)
    np.random.seed(seed)
    os.makedirs(merge_save_dir, exist_ok=True)

    circles = [((random.uniform(-4, 4), random.uniform(-4, 4)), random.uniform(1.0, 2.4)) for _ in range(n)]
    history = []  # record merges as dicts
    fig, ax = plt.subplots(figsize=(6, 6))
    step = 0

    # Save initial frame
    plot_circles_merge(ax, circles, step="merge init")
    plt.savefig(os.path.join(merge_save_dir, f"merge_step_{step:03d}_initial.png"), dpi=150)

    while len(circles) > 1:
        step += 1
        # find closest pair according to circle_distance
        best_i = best_j = None
        best_d = float('inf')
        for i in range(len(circles)):
            for j in range(i + 1, len(circles)):
                d = circle_distance(circles[i], circles[j])
                if d < best_d:
                    best_d, best_i, best_j = d, i, j

        # save frame highlighting pair about to merge
        plot_circles_merge(ax, circles, highlight=(best_i, best_j), step=step)
        plt.savefig(os.path.join(merge_save_dir, f"merge_step_{step:03d}_pair.png"), dpi=150)

        # form proxy and save merged frame
        c_i = circles[best_i]
        c_j = circles[best_j]
        proxy = make_combined_circle(c_i[0], c_i[1], c_j[0], c_j[1], deterministic=deterministic_proxy)
        plot_circles_merge(ax, circles, highlight=(best_i, best_j), proxy=proxy, step=f"{step} (merged)")
        plt.savefig(os.path.join(merge_save_dir, f"merge_step_{step:03d}_merged.png"), dpi=150)

        # record history (store the actual circles merged and proxy)
        history.append({'step': step, 'i': best_i, 'j': best_j, 'c_i': c_i, 'c_j': c_j, 'proxy': proxy})

        # remove originals and insert proxy
        for idx in sorted([best_i, best_j], reverse=True):
            circles.pop(idx)
        circles.append(proxy)

    # final frame (single root proxy)
    plot_circles_merge(ax, circles, step="merge final")
    plt.savefig(os.path.join(merge_save_dir, f"merge_step_final.png"), dpi=150)
    plt.close(fig)

    print(f"[merge] Saved {step*2 + 2} frames (including initial/final) to '{merge_save_dir}/'")
    # Return final proxy and history for the unmerge phase
    return circles[0], history

# ------------------------
# Unmerge / construction simulation
# ------------------------
def simulate_unmerge(history, root_proxy, seed=0, unmerge_save_dir="../unmerge_frames"):
    random.seed(seed)
    np.random.seed(seed)
    os.makedirs(unmerge_save_dir, exist_ok=True)

    tour = [{'pos': root_proxy[0], 'circles': [root_proxy]}]
    active_circles = [root_proxy]

    fig, ax = plt.subplots(figsize=(6, 6))
    frame_idx = 0

    def same_circle(a, b, eps=1e-9):
        return abs(a[0][0]-b[0][0]) < eps and \
               abs(a[0][1]-b[0][1]) < eps and \
               abs(a[1]-b[1]) < eps

    def dist(a, b):
        return math.hypot(a[0]-b[0], a[1]-b[1])

    # reverse merge history
    for entry in reversed(history):
        proxy = entry['proxy']
        c1, c2 = entry['c_i'], entry['c_j']

        removed_pts = []
        inserted_pts = []
        inserted_circles = []

        # --- remove proxy from tour ---
        for tp in tour:
            for circ in tp['circles']:
                if same_circle(circ, proxy):
                    tp['circles'].remove(circ)
                    removed_pts.append(tp['pos'])
                    break

        tour = [tp for tp in tour if tp['circles']]
        active_circles = [c for c in active_circles if not same_circle(c, proxy)]

        # --- insert generating circles ---
        for cgen in (c1, c2):
            assigned = False
            for tp in tour:
                if dist(tp['pos'], cgen[0]) <= cgen[1]:
                    tp['circles'].append(cgen)
                    assigned = True
                    break

            if not assigned:
                # approximate Alhazen insertion
                best_cost = float('inf')
                best_P = None
                best_idx = 0
                m = len(tour)
                O = cgen[0]

                def unit(v):
                    L = math.hypot(v[0], v[1])
                    return (v[0]/L, v[1]/L) if L > 1e-12 else (0.0, 0.0)

                for i in range(m):
                    A = tour[i]['pos']
                    B = tour[(i+1) % m]['pos']
                    uA = unit((A[0]-O[0], A[1]-O[1]))
                    uB = unit((B[0]-O[0], B[1]-O[1]))
                    bis = (uA[0]+uB[0], uA[1]+uB[1])
                    L = math.hypot(bis[0], bis[1])
                    if L < 1e-6:
                        continue
                    bis = (bis[0]/L, bis[1]/L)
                    P = (O[0] + 0.9*cgen[1]*bis[0],
                         O[1] + 0.9*cgen[1]*bis[1])
                    cost = dist(A,P) + dist(P,B) - dist(A,B)
                    if cost < best_cost:
                        best_cost = cost
                        best_P = P
                        best_idx = i+1

                # Fallback if no valid insertion point was found
                if best_P is None:
                    best_P = cgen[0]   # place tour point at circle center
                    best_idx = len(tour)

                tour.insert(best_idx, {'pos': best_P, 'circles': [cgen]})
                inserted_pts.append(best_P)


            active_circles.append(cgen)
            inserted_circles.append(cgen)

        # --- single frame for this unmerge ---
        plot_circles(
            ax,
            circles=active_circles,
            tour_pts=[tp['pos'] for tp in tour],
            removed_circles=[proxy],
            removed_pts=removed_pts,
            inserted_circles=inserted_circles,
            inserted_pts=inserted_pts,
            step=f"unmerge {entry['step']}"
        )

        plt.savefig(
            os.path.join(unmerge_save_dir,
                         f"unmerge_step_{frame_idx:03d}.png"),
            dpi=150
        )
        frame_idx += 1

    plt.close(fig)
    print(f"[unmerge] Saved {frame_idx} frames to '{unmerge_save_dir}/'")

# ------------------------
# Example usage (main)
# ------------------------
if __name__ == "__main__":
    # Configuration
    n = 8
    seed = 0
    merge_save_dir = "merge_save_dir"
    unmerge_save_dir = "unmerge_save_dir"

    # Run merge phase (saves merge frames and returns merge history)
    root_proxy, history = simulate_merge(n=n, seed=seed, merge_save_dir=merge_save_dir, deterministic_proxy=True)

    # Run unmerge phase using recorded history (saves unmerge frames)
    simulate_unmerge(history, root_proxy, seed=seed, unmerge_save_dir=unmerge_save_dir)

    print("Done. Merge frames in:", merge_save_dir)
    print("Done. Unmerge frames in:", unmerge_save_dir)
