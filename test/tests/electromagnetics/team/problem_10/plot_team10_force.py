#!/usr/bin/env python3
"""Plot the computed flux densities through S1, S2 and S3 against time.

Reads the postprocessor values written by team10.i and produces the figure
embedded in doc/content/problems/team_problem_10.md:

    team10_flux.png   the average flux density through each section against time

Fig. 2 of the problem statement defines the three sections: S1 across the centre
plate at z = 0, S2 across the top flange at x = 41.8 mm, and S3 across the web at
z = 0. B1 and B3 are the z component, B2 the x component. The flux runs down the
web, so team10.i reports B3 as negative; it is plotted here as -B3, in the
direction of the arrow in Fig. 2, which is the sense the benchmark tabulates.

Several runs can be overlaid for a convergence study, e.g. meshes with 2, 4 and
8 elements through the plate thickness:

    ./plot_team10_flux.py \\
        --computed OutputData/TEAM10CSV_2el.csv="2 elements" \\
        --computed OutputData/TEAM10CSV_4el.csv="4 elements" \\
        --computed OutputData/TEAM10CSV_8el.csv="8 elements"

The measurements are reported in the summary of results (Nakata, Takahashi and
Fujiwara, COMPEL 14, 1995). They are not reproduced here; supply them with
--measured as a CSV with columns section,time,B (section one of B1, B2, B3, time
in seconds, B in tesla), one row per point, so that each section can have its own
time samples as read off a published figure.

Usage:
    ./plot_team10_flux.py                               # computed only
    ./plot_team10_flux.py --measured team10_measured.csv
    ./plot_team10_flux.py --outdir /tmp
"""

import argparse
import csv
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]

DEFAULT_OUTDIR = REPO_ROOT / "doc" / "content" / "media"
DEFAULT_COMPUTED = HERE / "OutputData" / "TEAM10CSV.csv"

# Chart surface and ink, shared with the other TEAM figures.
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e6e5e1"
MEASURED = "#e8710a"
# One colour per computed run, in order; the first matches TEAM 20's.
COMPUTED = ["#2a78d6", "#1a9e77", "#7a5195", "#c43c5c"]

# Postprocessor name, the sign that turns it into the sense of Fig. 2's arrow,
# and the panel title.
SECTIONS = [
    ("B1", +1.0, "S1: centre plate, $B_z$"),
    ("B2", +1.0, "S2: top flange, $B_x$"),
    ("B3", -1.0, "S3: web, $-B_z$"),
]


def read_computed(path):
    """Return {section: (times, values)} from a team10.i postprocessor CSV."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        sys.exit(f"{path} has no rows")
    missing = [s for s, _, _ in SECTIONS if s not in rows[0]]
    if missing:
        sys.exit(f"{path} has no column(s) {', '.join(missing)}")
    t = [float(r["time"]) for r in rows]
    return {s: (t, [sign * float(r[s]) for r in rows]) for s, sign, _ in SECTIONS}


def read_measured(path):
    """Return {section: (times, values)} from a section,time,B CSV."""
    out = {s: ([], []) for s, _, _ in SECTIONS}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            s = r["section"].strip()
            if s not in out:
                sys.exit(f"{path}: unknown section {s!r}; expected B1, B2 or B3")
            out[s][0].append(float(r["time"]))
            out[s][1].append(float(r["B"]))
    for s, (t, b) in out.items():
        order = sorted(range(len(t)), key=t.__getitem__)
        out[s] = ([t[i] for i in order], [b[i] for i in order])
    return out


def interp(x, xs, ys):
    """Linear interpolation of (xs, ys) at x; None outside the sampled range."""
    if not xs or x < xs[0] or x > xs[-1]:
        return None
    for i in range(1, len(xs)):
        if x <= xs[i]:
            w = (x - xs[i - 1]) / (xs[i] - xs[i - 1]) if xs[i] > xs[i - 1] else 0.0
            return ys[i - 1] + w * (ys[i] - ys[i - 1])
    return ys[-1]


def _style(ax):
    """The axis styling used by the other figures in the documentation."""
    ax.set_facecolor(SURFACE)
    ax.set_axisbelow(True)
    ax.grid(True, axis="y", color=GRID, linewidth=1.0, solid_capstyle="butt")
    ax.grid(False, axis="x")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
        ax.spines[side].set_linewidth(1.0)
    ax.tick_params(which="both", colors=TEXT_SECONDARY, labelsize=9, length=0)


def plot(runs, measured, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(11.0, 4.0), dpi=200, sharex=True)
    fig.patch.set_facecolor(SURFACE)
    fig.suptitle("Average flux density through the measurement sections, TEAM Problem 10",
                 color=TEXT_PRIMARY, fontsize=12, x=0.008, ha="left", y=0.975)

    t_max = max(max(data[SECTIONS[0][0]][0]) for _, data in runs)
    for ax, (section, _, title) in zip(axes, SECTIONS):
        _style(ax)
        ax.set_title(title, color=TEXT_SECONDARY, fontsize=10, loc="left")
        if measured:
            t, b = measured[section]
            ax.plot(t, b, linestyle="none", marker="o", markersize=4.5, color=MEASURED,
                    markeredgecolor=SURFACE, markeredgewidth=1.0, zorder=3,
                    label="Measured")
        for (label, data), colour in zip(runs, COMPUTED):
            t, b = data[section]
            ax.plot(t, b, color=colour, linewidth=2.0, zorder=2, label=label)
        ax.set_xlim(0, t_max)
        ax.set_ylim(bottom=0)
        ax.set_xlabel("Time (s)", color=TEXT_SECONDARY, fontsize=10)
    axes[0].set_ylabel("Average flux density (T)", color=TEXT_SECONDARY, fontsize=10)

    legend = axes[0].legend(loc="lower right", frameon=False, fontsize=9)
    for text in legend.get_texts():
        text.set_color(TEXT_SECONDARY)

    fig.tight_layout(rect=(0, 0, 1, 0.93))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, facecolor=SURFACE)
    plt.close(fig)


def report(runs, measured):
    """Final values, and the RMS difference from the measurements where given."""
    for label, data in runs:
        print(label)
        for section, _, _ in SECTIONS:
            t, b = data[section]
            line = "  %s  at t = %.3f s: %.4f T" % (section, t[-1], b[-1])
            if measured and measured[section][0]:
                pairs = [(bm, interp(tm, t, b)) for tm, bm in zip(*measured[section])]
                diffs = [bc - bm for bm, bc in pairs if bc is not None]
                if diffs:
                    rms = math.sqrt(sum(d * d for d in diffs) / len(diffs))
                    line += "   RMS difference from measured %.4f T over %d points" % (
                        rms, len(diffs))
            print(line)


def parse_run(spec):
    """PATH or PATH=LABEL."""
    path, _, label = spec.partition("=")
    return Path(path), label or Path(path).stem


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--computed", action="append", metavar="PATH[=LABEL]",
                   help="a team10.i postprocessor CSV, optionally labelled; repeat to "
                        "overlay runs (default: %s)" % DEFAULT_COMPUTED.relative_to(HERE))
    p.add_argument("--measured", type=Path,
                   help="measured values as a section,time,B CSV (see above)")
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR,
                   help="directory the figure is written to (default: %(default)s)")
    args = p.parse_args()

    specs = args.computed or [f"{DEFAULT_COMPUTED}=Computed"]
    if len(specs) > len(COMPUTED):
        p.error(f"at most {len(COMPUTED)} computed runs can be overlaid")
    runs = []
    for spec in specs:
        path, label = parse_run(spec)
        if not path.exists():
            sys.exit(f"{path} not found; run team10.i first")
        runs.append((label, read_computed(path)))
    measured = read_measured(args.measured) if args.measured else None

    out = args.outdir / "team10_flux.png"
    plot(runs, measured, out)
    report(runs, measured)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
