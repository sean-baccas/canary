#!/usr/bin/env python3
"""Turn the TEAM Problem 10 B-H curve into expressions for the MOOSE input.

The benchmark specifies the steel by the (B, H) points of its Table 2, up to
B = 1.8 T and H = 9423 A/m, and above that by its Eq. (1),

    B = mu0 H + (a H^2 + b H + c)    1.8 <= B <= 2.22 T
    B = mu0 H + Ms                   B >= 2.22 T

with a = -2.381e-10, b = 2.327e-5, c = 1.590 and Ms = 2.16 T.  Below 1.8 T a
monotone cubic (Fritsch-Carlson) is put through the tabulated points, as for
TEAM Problem 20.  Above it, Eq. (1) is inverted exactly: the quadratic branch
has a closed-form H(B), and the saturation line is linear.

The solver needs the reluctivity nu = H/B and, for the Jacobian, B dnu/dB.  Both
follow from H(B) and its derivative,

    nu = H / B,        B dnu/dB = H'(B) - H/B,

so this script emits H(B) and H'(B) and the input file forms the rest.  Run it
with no arguments to print the two expressions, or --check to report how well
they reproduce the specification.

Two liberties are taken with the specification, both small, and --check
reports them:

* The quadratic gives B = 1.79997 T at the last tabulated field, H = 9423 A/m,
  not 1.8 T.  The last spline knot is moved to the quadratic's value so that the
  curve is continuous there, and the spline's end slope is clamped to the
  quadratic's so that it is continuously differentiable too.

* The magnetisation a H^2 + b H + c peaks at H* = -b / 2a = 48866 A/m, where it
  is 2.1586 T rather than the quoted Ms = 2.16 T, and B there is 2.2200 T.  Ms is
  taken as that peak, and the switch to the saturation line is made at H*, so the
  join is continuous with slope 1/mu0 on both sides.
"""

import argparse
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]
DEFAULT_OUTDIR = REPO_ROOT / "doc" / "content" / "media"

# Chart surface and ink, shared with the other TEAM figures.
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e6e5e1"
CURVE = "#2a78d6"

MU0 = 4.0e-7 * math.pi

# B (T), H (A/m) from Table 2 of the problem statement.  The table prints
# (0.8, 289) at the foot of one column and the head of the next; it is listed
# once here.
TABLE = [
    (0.0, 0.0), (0.0025, 16.0), (0.005, 30.0), (0.0125, 54.0), (0.025, 93.0),
    (0.05, 143.0), (0.1, 191.0), (0.2, 210.0), (0.3, 222.0), (0.4, 233.0),
    (0.5, 247.0), (0.6, 258.0), (0.7, 272.0), (0.8, 289.0), (0.9, 313.0),
    (1.00, 342.0), (1.10, 377.0), (1.20, 433.0), (1.30, 509.0), (1.40, 648.0),
    (1.50, 933.0), (1.55, 1228.0), (1.60, 1934.0), (1.65, 2913.0),
    (1.70, 4993.0), (1.75, 7189.0), (1.80, 9423.0),
]

# Eq. (1) of the problem statement.
EQ_A, EQ_B, EQ_C = -2.381e-10, 2.327e-5, 1.590
MS_QUOTED = 2.16
B_QUOTED_SWITCH = 2.22

BETA = EQ_B + MU0  # dB/dH of the quadratic branch at H = 0


def b_quadratic(h):
    """B(H) on the quadratic branch of Eq. (1)."""
    return MU0 * h + EQ_A * h * h + EQ_B * h + EQ_C


def h_quadratic(b):
    """The rising root of b = a H^2 + beta H + c, in a form free of cancellation."""
    return 2.0 * (b - EQ_C) / (BETA + math.sqrt(BETA * BETA + 4.0 * EQ_A * (b - EQ_C)))


def dh_quadratic(b):
    """dH/dB on the quadratic branch: 1 / (2 a H + beta) = 1 / sqrt(discriminant)."""
    return 1.0 / math.sqrt(BETA * BETA + 4.0 * EQ_A * (b - EQ_C))


H_LAST = TABLE[-1][1]
B_LAST = b_quadratic(H_LAST)  # 1.79997 T: the last knot, moved onto Eq. (1)
H_STAR = -EQ_B / (2.0 * EQ_A)  # 48866 A/m: where the magnetisation peaks
M_SAT = EQ_C - EQ_B * EQ_B / (4.0 * EQ_A)  # 2.1586 T: the peak magnetisation
B_SWITCH = MU0 * H_STAR + M_SAT  # 2.2200 T: where the saturation line takes over

KNOTS = TABLE[:-1] + [(B_LAST, H_LAST)]


def slopes(bs, hs):
    """Fritsch-Carlson slopes, clamped at the right end onto the quadratic's slope."""
    n = len(bs)
    h = [bs[i + 1] - bs[i] for i in range(n - 1)]
    d = [(hs[i + 1] - hs[i]) / h[i] for i in range(n - 1)]
    m = [0.0] * n
    for i in range(1, n - 1):
        if d[i - 1] * d[i] > 0.0:
            w1, w2 = 2 * h[i] + h[i - 1], h[i] + 2 * h[i - 1]
            m[i] = (w1 + w2) / (w1 / d[i - 1] + w2 / d[i])
    m[0] = d[0]
    m[-1] = min(dh_quadratic(bs[-1]), 3.0 * d[-1])
    return m


def _n(v):
    return "%.7g" % v


def _segments():
    bs = [b for b, _ in KNOTS]
    hs = [h for _, h in KNOTS]
    m = slopes(bs, hs)
    for i in range(len(bs) - 1):
        w = bs[i + 1] - bs[i]
        c2 = (3 * (hs[i + 1] - hs[i]) / w - 2 * m[i] - m[i + 1]) / w
        c3 = (m[i] + m[i + 1] - 2 * (hs[i + 1] - hs[i]) / w) / (w * w)
        yield bs[i], bs[i + 1], hs[i], m[i], c2, c3


def expressions():
    """Return the H(B) and dH/dB expressions, in terms of a symbol named B."""
    h_segs, d_segs = [], []
    for b0, b1, h0, m0, c2, c3 in _segments():
        t = "B" if b0 == 0.0 else "(B-%s)" % _n(b0)
        h_segs.append((b1, "%s+%s*(%s+%s*(%s+%s*%s))"
                       % (_n(h0), t, _n(m0), t, _n(c2), t, _n(c3))))
        d_segs.append((b1, "%s+%s*(%s+%s*%s)" % (_n(m0), t, _n(2 * c2), t, _n(3 * c3))))

    # Eq. (1): the quadratic branch up to B_SWITCH, then the saturation line.
    disc = "sqrt(%s+%s*(B-%s))" % (_n(BETA * BETA), _n(4.0 * EQ_A), _n(EQ_C))
    h_tail = "if(B<%s,2*(B-%s)/(%s+%s),(B-%s)/%s)" % (
        _n(B_SWITCH), _n(EQ_C), _n(BETA), disc, _n(M_SAT), _n(MU0))
    d_tail = "if(B<%s,1/%s,%s)" % (_n(B_SWITCH), disc, _n(1.0 / MU0))

    def wrap(segs, tail):
        e = tail
        for bmax, seg in reversed(segs):
            e = "if(B<%s,%s,%s)" % (_n(bmax), seg, e)
        return e

    return wrap(h_segs, h_tail), wrap(d_segs, d_tail)


# --- evaluation, used only by --check and --plot ----------------------------
def h_of_b(b):
    for b0, b1, h0, m0, c2, c3 in _segments():
        if b < b1:
            t = b - b0
            return h0 + t * (m0 + t * (c2 + t * c3))
    if b < B_SWITCH:
        return h_quadratic(b)
    return (b - M_SAT) / MU0


def check():
    print("Eq.(1) gives B = %.5f T at the last tabulated field, H = %.0f A/m; the "
          "table says %.2f" % (B_LAST, H_LAST, TABLE[-1][0]))
    print("so the last spline knot is moved to B = %.5f T\n" % B_LAST)
    print("the magnetisation of Eq.(1) peaks at H = %.0f A/m with M = %.4f T, B = %.4f T;"
          % (H_STAR, M_SAT, B_SWITCH))
    print("the statement quotes Ms = %.2f T and B = %.2f T, so Ms is taken as %.4f T "
          "and the switch made at %.4f T\n" % (MS_QUOTED, B_QUOTED_SWITCH, M_SAT, B_SWITCH))

    worst = max(abs(h_of_b(b) - h) for b, h in KNOTS[:-1])
    print("tabulated points reproduced to %.3g A/m" % worst)

    n, prev, ok = 400000, -1.0, True
    for i in range(n + 1):
        v = h_of_b(2.6 * i / n)
        if v < prev:
            ok = False
        prev = v
    print("H(B) monotone increasing on [0, 2.6]: %s" % ok)

    eps = 1e-7
    for name, bj in (("spline/quadratic", B_LAST), ("quadratic/saturation", B_SWITCH)):
        left = (h_of_b(bj - eps) - h_of_b(bj - 2 * eps)) / eps
        right = (h_of_b(bj + 2 * eps) - h_of_b(bj + eps)) / eps
        jump = h_of_b(bj + 1e-12) - h_of_b(bj - 1e-12)
        print("%s join at %.5f T: jump %.2g A/m, slopes %.0f and %.0f A/m/T"
              % (name, bj, jump, left, right))
    print("(1/mu0 = %.0f A/m/T)" % (1.0 / MU0))

    nu = [(h_of_b(b) / b, b) for b in [0.001 * k for k in range(1, 2400)]]
    lo = min(nu)
    print("\nreluctivity nu = H/B runs from %.0f as B -> 0, down to %.0f at B = %.2f T,"
          % (h_of_b(1e-6) / 1e-6, lo[0], lo[1]))
    print("and up to %.0f at B = %.2f T" % (h_of_b(2.4) / 2.4, 2.4))


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


def plot(output):
    """B against H on a log field axis, over the whole range through to saturation."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(figsize=(7.6, 4.6), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    _style(ax)
    fig.suptitle("B-H curve of the TEAM Problem 10 steel", color=TEXT_PRIMARY,
                 fontsize=12, x=0.012, ha="left", y=0.975)
    ax.set_xscale("log")

    n, b_max = 6000, 2.4
    bs = [b_max * i / n for i in range(1, n + 1)]
    spline = [(h_of_b(b), b) for b in bs if b <= B_LAST]
    tail = [(h_of_b(b), b) for b in [B_LAST] + [b for b in bs if b > B_LAST]]
    ax.plot([h for h, _ in spline], [b for _, b in spline], color=CURVE, linewidth=2.0)
    ax.plot([h for h, _ in tail], [b for _, b in tail], color=CURVE, linewidth=2.0,
            linestyle=(0, (4, 2)))
    ax.plot([h for _, h in TABLE[1:]], [b for b, _ in TABLE[1:]], linestyle="none",
            marker="o", markersize=5.0, color=CURVE, markeredgecolor="none")
    ax.set_xlim(10, 1.0e6)
    ax.set_ylim(0, b_max)
    ax.set_xlabel("Magnetic field strength, H (A/m)", color=TEXT_SECONDARY, fontsize=10)
    ax.set_ylabel("Flux density, B (T)", color=TEXT_SECONDARY, fontsize=10)
    fig.text(0.988, 0.915, "Eq. (1) above B = 1.8 T (dashed), $M_s$ = %.4f T" % M_SAT,
             color=TEXT_SECONDARY, fontsize=8.5, ha="right", va="top")
    handles = [
        Line2D([], [], linestyle="none", marker="o", markersize=5.0, color=CURVE,
               markeredgecolor="none", label="Table 2 of the problem statement"),
        Line2D([], [], color=CURVE, linewidth=2.0, label="Monotone cubic spline"),
    ]
    legend = ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=9)
    for text in legend.get_texts():
        text.set_color(TEXT_SECONDARY)

    fig.tight_layout(rect=(0, 0, 1, 0.9))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, facecolor=SURFACE)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--check", action="store_true",
                   help="report how faithfully the expressions follow the specification")
    p.add_argument("--plot", action="store_true",
                   help="write a figure of the B-H curve instead of the expressions")
    p.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR,
                   help="directory the figure is written to (default: %(default)s)")
    args = p.parse_args()
    if args.check:
        check()
        return
    if args.plot:
        out = args.outdir / "team10_bh_curve.png"
        plot(out)
        print(f"wrote {out}")
        return
    h, d = expressions()
    print("# H(B), A/m -- paste as the expression of the h_of_b function")
    print(h)
    print()
    print("# dH/dB, A/m/T -- paste as the expression of the dh_db function")
    print(d)


if __name__ == "__main__":
    main()
