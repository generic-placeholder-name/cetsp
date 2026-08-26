#!/usr/bin/env python3
"""
grad.py

Projected gradient-descent (Python translation of the provided C++ logic).
Saves plot to ../grad.png

Usage: python grad.py
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

EPS = 1e-12

# ----------------------------
# Geometry / helpers
# ----------------------------

def angle_at_P(P, A, B):
    u = A - P
    v = B - P
    nu = np.linalg.norm(u)
    nv = np.linalg.norm(v)
    if nu < EPS or nv < EPS:
        return 0.0
    cos = np.dot(u, v) / (nu * nv)
    cos = np.clip(cos, -1.0, 1.0)
    return math.acos(cos)

def in_all_disks(pt, circles, tol=1e-12):
    for (cx, cy, r) in circles:
        if np.hypot(pt[0]-cx, pt[1]-cy) > r + tol:
            return False
    return True

# ----------------------------
# Projected gradient step (faithful to the C++ logic)
# ----------------------------

def gradient_and_max_step(p, a, b, circles):
    """
    Given p (2-array), a, b (2-arrays), and circles list of (center(2), r),
    compute gradient direction (normalized sum of unit vectors to a and b),
    and compute maxStep before leaving any circle using quadratic solve as in C++.
    Returns (d (unit vec) or None, maxStep float).
    """
    x, y = p[0], p[1]

    # accumulate gradient components (sum of unit directions toward a and b)
    gx = 0.0
    gy = 0.0

    # toward A
    dx = a[0] - x
    dy = a[1] - y
    d = math.hypot(dx, dy)
    if d > EPS:
        gx += dx / d
        gy += dy / d

    # toward B
    dx = b[0] - x
    dy = b[1] - y
    d = math.hypot(dx, dy)
    if d > EPS:
        gx += dx / d
        gy += dy / d

    norm = math.hypot(gx, gy)
    if norm < EPS:
        # degenerate gradient
        return None, 0.0

    # normalize to unit step direction
    gx /= norm
    gy /= norm

    maxStep = float('inf')
    # for each circle solve a2*t^2 + b2*t + c2 = 0
    for (cx, cy, r) in circles:
        dx = x - cx
        dy = y - cy
        a2 = gx*gx + gy*gy      # should be 1.0 due to normalization, but keep general
        b2 = 2.0 * (dx*gx + dy*gy)
        c2 = dx*dx + dy*dy - r*r

        # discriminant
        disc = b2*b2 - 4.0*a2*c2
        if disc < 0.0:
            # no real roots -> ray misses circle boundary (numerically possible)
            # In C++ they clamp disc to 0 via max(..., 0.)
            # here we skip negative
            continue
        disc = max(disc, 0.0)
        sqrt_disc = math.sqrt(disc)
        # two roots
        t1 = (-b2 + sqrt_disc) / (2.0 * a2)
        t2 = (-b2 - sqrt_disc) / (2.0 * a2)
        tmax = max(t1, t2, 0.0)   # faithful to provided C++: take max(t1,t2,0)
        if tmax < maxStep:
            maxStep = tmax

    return np.array([gx, gy], dtype=float), maxStep

def descent_max(circles, p0, a, b, steps=10):
    p = p0.astype(float).copy()
    path = [p.copy()]
    for _ in range(steps):
        d, t = gradient_and_max_step(p, a, b, circles)
        if d is None or t <= 0.0 or not np.isfinite(t):
            break
        # advance
        p = p + t * d
        path.append(p.copy())
    return np.array(path)

# ----------------------------
# Example config generation
# ----------------------------

def generate_example():
    # deterministic overlapping circles
    circles = [
        (0.0, 0.0, 1.2),
        (0.8, 0.2, 1.05),
        (0.35, 0.9, 0.95),
    ]
    # pick P inside intersection
    P = np.array([0.42, 0.35], dtype=float)
    assert in_all_disks(P, circles), "Generated P not inside all circles"

    # choose A, B such that AB doesn't cross intersection (sample check)
    A = np.array([2.2, 1.6], dtype=float)
    B = np.array([-1.8, 1.6], dtype=float)
    # crude check: sample segment AB
    samples = np.linspace(0.0, 1.0, 401)[:,None] * (B - A) + A
    intersects = any(in_all_disks(pt, circles) for pt in samples)
    if intersects:
        print("Warning: AB intersects intersection in the generated config.")
    return circles, A, B, P

# ----------------------------
# Plotting
# ----------------------------

def plot_scene(circles, A, B, path, filename="../grad.png"):
    fig, ax = plt.subplots(figsize=(9,7))
    ax.set_aspect('equal', 'box')

    # draw circles
    for (cx, cy, r) in circles:
        circ = Circle((cx, cy), r, fill=False, lw=1.4, alpha=0.9)
        ax.add_patch(circ)

    # fill intersection (approx)
    xs = [cx - r for (cx,_,r) in circles] + [cx + r for (cx,_,r) in circles]
    ys = [cy - r for (_,cy,r) in circles] + [cy + r for (_,cy,r) in circles]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    xmin, xmax = min(xmin, A[0], B[0]) - 0.5, max(xmax, A[0], B[0]) + 0.5
    ymin, ymax = min(ymin, A[1], B[1]) - 0.5, max(ymax, A[1], B[1]) + 0.5
    nx, ny = 300, 300
    gx = np.linspace(xmin, xmax, nx)
    gy = np.linspace(ymin, ymax, ny)
    X, Y = np.meshgrid(gx, gy)
    mask = np.ones_like(X, dtype=bool)
    for (cx, cy, r) in circles:
        mask &= ((X-cx)**2 + (Y-cy)**2 <= r*r + 1e-12)
    ax.contourf(X, Y, mask.astype(float), levels=[0.5, 1.5], alpha=0.22, colors=['#bfefff'])

    # extract P and P'
    path = np.asarray(path)
    P = path[0]
    Pp = path[-1]

    # --- Segments ---
    # Original: A–P–B
    ax.plot([A[0], P[0]], [A[1], P[1]], 'k--', lw=1.2, label='A–P–B (orig)')
    ax.plot([P[0], B[0]], [P[1], B[1]], 'k--', lw=1.2)
    # New: A–P'–B
    ax.plot([A[0], Pp[0]], [A[1], Pp[1]], 'm-', lw=1.4, label="A–P'–B (new)")
    ax.plot([Pp[0], B[0]], [Pp[1], B[1]], 'm-', lw=1.4)

    # --- Points ---
    ax.plot(A[0], A[1], 'ro', label='A')
    ax.plot(B[0], B[1], 'go', label='B')
    ax.plot(P[0], P[1], 'ko', label='P')
    ax.plot(Pp[0], Pp[1], 'mo', label="P'")

    # --- Labels inside graph ---
    ax.text(A[0]+0.05, A[1]+0.05, 'A', color='r', fontsize=10, fontweight='bold')
    ax.text(B[0]+0.05, B[1]+0.05, 'B', color='g', fontsize=10, fontweight='bold')
    ax.text(P[0]+0.05, P[1], 'P', color='k', fontsize=10, fontweight='bold')
    ax.text(Pp[0]+0.05, Pp[1], "P'", color='m', fontsize=10, fontweight='bold')

    # --- Optional path (if multiple steps) ---
    if len(path) > 1:
        ax.plot(path[:,0], path[:,1], 'k:', lw=0.8, alpha=0.6)

    ax.legend(loc='upper right')
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_title("Projected Gradient Descent Reoptimization")
    plt.tight_layout()
    plt.savefig(filename, dpi=200)
    plt.close(fig)
    print(f"Saved plot to {filename}")

# ----------------------------
# Main
# ----------------------------

def main():
    circles, A, B, P = generate_example()
    print("Circles:", circles)
    print("A:", A, "B:", B, "P start:", P)
    print("Initial angle ∠APB (deg):", math.degrees(angle_at_P(P, A, B)))

    path = descent_max(circles, P, A, B, steps=20)
    P_end = path[-1]
    print("P end:", P_end)
    print("Final angle ∠APB (deg):", math.degrees(angle_at_P(P_end, A, B)))

    plot_scene(circles, A, B, path, filename="../grad.png")

if __name__ == "__main__":
    main()
