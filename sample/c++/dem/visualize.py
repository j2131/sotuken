#!/usr/bin/env python3
# /// script
# dependencies = [
#   "numpy",
#   "matplotlib",
# ]
# ///
import argparse
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import animation

# Cylinder size (meters)
CYL_RADIUS = 0.005   # 0.5 cm
CYL_LENGTH = 0.05    # 5 cm
PLOT_MARGIN = 0.002  # 2 mm margin for visibility

def load_dat(path):
    with open(path, "r", encoding="utf-8") as f:
        lines = f.read().strip().splitlines()
    if len(lines) < 2:
        return None, None
    time = float(lines[0].strip())
    nbody = int(lines[1].strip())
    data = []
    for line in lines[2:2 + nbody]:
        cols = line.strip().split()
        if len(cols) < 15:
            continue
        # id(0) mass(1) rad(2) x(3) y(4) z(5) ... kind(15)
        # Check if kind column exists, else default to 0
        p_kind = int(cols[15]) if len(cols) > 15 else 0
        data.append([float(cols[2]), float(cols[3]), float(cols[4]), float(cols[5]), p_kind])
    if not data:
        return time, None
    arr = np.array(data)
    rad = arr[:, 0]
    pos = arr[:, 1:4]
    kind = arr[:, 4].astype(int)
    return time, (pos, rad, kind)

def main():
    parser = argparse.ArgumentParser(description="Visualize DEM results (*.dat)")
    parser.add_argument("result_dir", nargs="?", default="result", help="result directory")
    parser.add_argument("--fps", type=int, default=30, help="animation fps")
    parser.add_argument("--index", type=int, default=None, help="show single frame index")
    parser.add_argument("--save", type=str, default=None, help="save animation to file (e.g. out.mp4 or out.gif)")
    parser.add_argument("--mode", choices=["3d", "2d"], default="3d", help="plot mode")
    parser.add_argument("--view", choices=["x", "y", "z"], default="z", help="view axis for 2d (look along axis)")
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.result_dir, "*.dat")))
    if not files:
        raise SystemExit(f"No .dat files found in {args.result_dir}")

    half_len = CYL_LENGTH * 0.5 + PLOT_MARGIN
    rad_lim = CYL_RADIUS + PLOT_MARGIN

    if args.index is not None:
        idx = max(0, min(args.index, len(files) - 1))
        t, payload = load_dat(files[idx])
        if payload is None:
            raise SystemExit("No particle data")
        pos, rad, kind = payload
        colors = ["red", "green", "blue"]
        c_list = [colors[k] for k in kind]

        if args.mode == "2d":
            fig = plt.figure(figsize=(6, 6))
            ax = fig.add_subplot(111)
            if args.view == "x":
                xy = pos[:, 1:3]
                ax.set_xlabel("y")
                ax.set_ylabel("z")
            elif args.view == "y":
                xy = pos[:, [0, 2]]
                ax.set_xlabel("x")
                ax.set_ylabel("z")
            else:
                xy = pos[:, 0:2]
                ax.set_xlabel("x")
                ax.set_ylabel("y")
            for i in range(len(rad)):
                circ = Circle((xy[i, 0], xy[i, 1]), rad[i], color=c_list[i], alpha=0.7)
                ax.add_patch(circ)
            ax.set_aspect("equal", adjustable="box")
            ax.set_title(f"t={t:.4f}")
            if args.view == "x":
                ax.set_xlim(-half_len, half_len)
                ax.set_ylim(-rad_lim, rad_lim)
            elif args.view == "y":
                ax.set_xlim(-rad_lim, rad_lim)
                ax.set_ylim(-rad_lim, rad_lim)
            else:
                ax.set_xlim(-rad_lim, rad_lim)
                ax.set_ylim(-half_len, half_len)
        else:
            fig = plt.figure(figsize=(8, 8))
            ax = fig.add_subplot(111, projection="3d")
            # Scatter particles with their actual radius (scaling factor 50000 for visibility of small particles)
            ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], s=(rad * 50000) ** 2, c=c_list, alpha=0.7)
            ax.set_title(f"t={t:.4f}")
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_zlabel("z")
            # Set bounds to see the whole cylinder area
            ax.set_xlim(-rad_lim, rad_lim)
            ax.set_ylim(-half_len, half_len)
            ax.set_zlim(-rad_lim, rad_lim)

        if args.save:
            plt.savefig(args.save)
            print(f"Saved image to {args.save}")
        else:
            plt.show()
        return

    frames = []
    times = []
    for path in files:
        t, payload = load_dat(path)
        if payload is None:
            continue
        pos, rad, kind = payload
        frames.append((pos, rad, kind))
        times.append(t)

    colors = ["red", "green", "blue"]
    if args.mode == "2d":
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)
        ax.set_aspect("equal", adjustable="box")
        if args.view == "x":
            ax.set_xlim(-half_len, half_len)
            ax.set_ylim(-rad_lim, rad_lim)
        elif args.view == "y":
            ax.set_xlim(-rad_lim, rad_lim)
            ax.set_ylim(-rad_lim, rad_lim)
        else:
            ax.set_xlim(-rad_lim, rad_lim)
            ax.set_ylim(-half_len, half_len)
        if args.view == "x":
            ax.set_xlabel("y")
            ax.set_ylabel("z")
        elif args.view == "y":
            ax.set_xlabel("x")
            ax.set_ylabel("z")
        else:
            ax.set_xlabel("x")
            ax.set_ylabel("y")
        patches = []
        for _ in range(len(frames[0][1])):
            circ = Circle((0.0, 0.0), 0.0, alpha=0.7)
            ax.add_patch(circ)
            patches.append(circ)
    else:
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection="3d")
        scat = ax.scatter([], [], [], s=[], alpha=0.7)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        # Set bounds to see the whole cylinder area
        ax.set_xlim(-rad_lim, rad_lim)
        ax.set_ylim(-half_len, half_len)
        ax.set_zlim(-rad_lim, rad_lim)

    def init():
        if args.mode == "2d":
            for p in patches:
                p.set_radius(0.0)
            return patches
        else:
            scat._offsets3d = ([], [], [])
            scat.set_sizes([])
            return scat,

    def update(i):
        pos, rad, kind = frames[i]
        if args.mode == "2d":
            if args.view == "x":
                xy = pos[:, 1:3]
            elif args.view == "y":
                xy = pos[:, [0, 2]]
            else:
                xy = pos[:, 0:2]
            for j, p in enumerate(patches):
                p.center = (xy[j, 0], xy[j, 1])
                p.set_radius(rad[j])
                p.set_color(colors[kind[j]])
            ax.set_title(f"t={times[i]:.4f}")
            return patches
        else:
            scat._offsets3d = (pos[:, 0], pos[:, 1], pos[:, 2])
            scat.set_sizes((rad * 50000) ** 2)
            scat.set_color([colors[k] for k in kind])
            ax.set_title(f"t={times[i]:.4f}")
            return scat,

    ani = animation.FuncAnimation(
        fig, update, frames=len(frames), init_func=init, interval=1000 / args.fps, blit=False
    )

    if args.save:
        print(f"Saving animation to {args.save}...")
        if args.save.endswith(".gif"):
            ani.save(args.save, writer='pillow', fps=args.fps)
        else:
            ani.save(args.save, writer='ffmpeg', fps=args.fps)
        print("Done.")
    else:
        plt.show()

if __name__ == "__main__":
    main()
