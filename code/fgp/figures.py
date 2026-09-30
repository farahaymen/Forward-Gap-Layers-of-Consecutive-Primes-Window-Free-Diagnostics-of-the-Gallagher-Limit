"""
fgp.figures -- every figure in the paper, in one style.

Layout policy (the reason the earlier version produced unreadable figures):

  * Panel labels are drawn OUTSIDE the axes frame, as a left-aligned title.
    Drawing them inside at the top-left corner put them exactly where a legend
    placed at "upper left" lands, so the two overlapped.
  * Legends are drawn OUTSIDE the data area, as figure-level legends beneath
    the panels ("outside lower center", which constrained layout reserves room
    for). No in-axes legend can then cover a curve, whatever the data does.
  * A family of twenty curves indexed by scale is a sequential encoding, so it
    gets a colourbar rather than a legend listing three of the twenty.
  * Axis limits are computed from the data with a fixed fractional pad, so a
    rerun on a different range cannot clip a series.

Figure widths match the text width of the manuscript (A4, 2.6 cm margins ->
6.22 in), so \\includegraphics[width=\\linewidth] scales them 1:1 and the font
sizes below are the sizes that appear on the page.

Colours are the Okabe-Ito colour-blind-safe set; the categorical subsets used
here pass adjacent-pair CVD separation, and every series additionally differs
by marker or line style, so identity is never carried by colour alone.

All figures are written as PDF (vector, for the manuscript) and PNG (preview).
"""

import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from .spatial import R_model

# ------------------------------------------------------------------ style
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Times New Roman", "Times"],
    "mathtext.fontset": "dejavuserif",
    "font.size": 8.0,
    "axes.labelsize": 8.0,
    "axes.titlesize": 8.0,
    "legend.fontsize": 7.0,
    "xtick.labelsize": 7.0,
    "ytick.labelsize": 7.0,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.45,
    "ytick.minor.width": 0.45,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "lines.linewidth": 1.0,
    "legend.frameon": False,
    "legend.handlelength": 2.2,
    "legend.columnspacing": 1.4,
    "legend.handletextpad": 0.5,
    "figure.dpi": 150,
    "savefig.dpi": 400,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# Okabe-Ito colour-blind-safe palette
C = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73",
     "red": "#D55E00", "purple": "#CC79A7", "sky": "#56B4E9",
     "yellow": "#F0E442", "black": "#000000", "grey": "#6e6e6e"}

# Text width of the manuscript (A4 paper, 2.6 cm margins) in inches.
FULL = 6.22
HALF = 3.05
SEQ_CMAP = "viridis"          # perceptually uniform, monotone in lightness


def _save(fig, outdir, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"),
                    bbox_inches="tight", pad_inches=0.015)
    plt.close(fig)


def _panel(ax, tag, text=""):
    """Panel label above the frame, so nothing inside the axes can collide."""
    ax.set_title(f"{tag} {text}".rstrip(), loc="left", fontsize=8.0, pad=3.5)


def _pad_ylim(ax, values, frac=0.08, lo=None, hi=None):
    """Set y limits from the data with a fractional pad on each side."""
    v = np.asarray([x for x in np.ravel(values) if np.isfinite(x)], dtype=float)
    if v.size == 0:
        return
    a, b = float(v.min()), float(v.max())
    span = (b - a) or max(abs(b), 1.0) * 0.1
    ax.set_ylim(a - frac * span if lo is None else lo,
                b + frac * span if hi is None else hi)


def _below(fig, handles, labels, ncol, fontsize=7.0):
    """One figure-level legend under the panels.

    No bbox_to_anchor: constrained layout reserves the strip itself when the
    location is one of the "outside" forms, and an explicit anchor overrides
    that reservation and drops the legend back on top of the tick labels.
    """
    return fig.legend(handles, labels, loc="outside lower center", ncol=ncol,
                      fontsize=fontsize, borderaxespad=0.4)


# ------------------------------------------------------------ Figure 1
def fig_cfit(prim: pd.DataFrame, outdir, quoted=None):
    """C_fit under the primary specification, the envelope over all
    specifications, and the Wolf benchmark, against X.

    `quoted` is an optional (X, value, label) triple for a value reported in the
    literature. It is omitted unless the caller supplies one, so the figure never
    asserts an attribution the manuscript has not verified.
    """
    fig, ax = plt.subplots(figsize=(FULL, 2.45), layout="constrained")
    X = prim.X.to_numpy()

    ax.fill_between(X, prim.env_min, prim.env_max, color=C["sky"], alpha=0.30,
                    lw=0, label="envelope of all specifications", zorder=1)
    ax.axhline(1, color=C["black"], lw=0.7, ls=":", zorder=2,
               label="Gallagher limit")
    ax.plot(X, prim.C_W, "s--", ms=2.8, lw=0.9, color=C["orange"], zorder=3,
            label=r"$C_W(X)=\pi(X)\log X/X$")
    ax.errorbar(X, prim.C_fit, yerr=2 * prim.SE, fmt="o-", ms=3, lw=1.0,
                capsize=1.5, color=C["blue"], zorder=4,
                label=r"$C_{\mathrm{fit}}(X)$, primary spec. ($\pm 2$ SE)")
    if quoted is not None:
        qx, qy, qlab = quoted
        ax.plot([qx], [qy], marker="*", ms=9, color=C["red"], ls="none",
                zorder=5, label=qlab)

    ax.set_xscale("log")
    ax.set_xlabel(r"$X = p_N$")
    ax.set_ylabel("rescaled decay rate")
    vals = [prim.env_min, prim.env_max, prim.C_W, prim.C_fit, [1.0]]
    if quoted is not None:
        vals.append([quoted[1]])
    _pad_ylim(ax, np.concatenate([np.asarray(v, dtype=float).ravel() for v in vals]))
    h, l = ax.get_legend_handles_labels()
    _below(fig, h, l, ncol=4 if quoted is None else 5)
    _save(fig, outdir, "fig1_cfit")


# ------------------------------------------------------------ Figure 2
def fig_specgrid(grid: pd.DataFrame, outdir):
    """Every specification across every scale: colour = window, style = estimator."""
    fig, axes = plt.subplots(1, 2, figsize=(FULL, 2.55), sharey=True,
                             layout="constrained")
    styles = {"POIS": "-", "WLS": "--", "OLS": ":"}
    # The window width c is an ordered variable, not a set of categories, so it
    # is encoded by a sequential ramp (narrow = light, wide = dark) rather than
    # by unrelated hues. That also makes the figure's message legible directly:
    # widening the window raises the fitted rate.
    cs = sorted(grid.c.unique())
    ramp = plt.get_cmap(SEQ_CMAP)(np.linspace(0.72, 0.02, len(cs)))
    cols = {c: ramp[i] for i, c in enumerate(cs)}

    for ax, gmin in zip(axes, sorted(grid.g_min.unique())):
        sub = grid[grid.g_min == gmin]
        for (c, est), grp in sub.groupby(["c", "estimator"]):
            grp = grp.sort_values("X")
            ax.plot(grp.X, grp.C_fit, styles[est], color=cols[c], lw=0.9)
        ax.axhline(1, color=C["black"], lw=0.7, ls=":", zorder=0)
        ax.set_xscale("log")
        ax.set_xlabel(r"$X = p_N$")

    _panel(axes[0], "(a)", r"$g \geq 2$")
    _panel(axes[1], "(b)", r"$g \geq 4$")
    axes[0].set_ylabel(r"$C_{\mathrm{fit}}(X)$")
    _pad_ylim(axes[0], np.r_[grid.C_fit.to_numpy(), 1.0])

    win = [Line2D([], [], color=cols[c], lw=1.6,
                  label=rf"$g \leq {c}\log X$") for c in cs]
    est = [Line2D([], [], color=C["grey"], lw=1.4, ls=styles[e], label=lab)
           for e, lab in (("POIS", "Poisson GLM"), ("WLS", "weighted LS"),
                          ("OLS", "unweighted LS"))]
    _below(fig, win + est, [h.get_label() for h in win + est], ncol=5)
    _save(fig, outdir, "fig2_specgrid")


# ------------------------------------------------------------ Figure 3
def fig_local(loc: pd.DataFrame, pl: pd.DataFrame, outdir):
    """The window-free local diagnostics against scale, with their power laws.

    The fitted exponent goes in the panel title rather than an in-axes
    annotation, so it cannot land on the data.
    """
    fig, axes = plt.subplots(1, 3, figsize=(FULL, 2.35), layout="constrained")
    L = loc.logX_eff.to_numpy()
    Lf = np.linspace(L.min() - 0.25, L.max() + 0.25, 200)
    plr = pl.set_index("statistic")

    panels = [("E2", "SE_E2", 1.0, r"$E_2 = m_2/(2m_1^2)$", "(a)"),
              ("E3", "SE_E3", 1.0, r"$E_3 = m_3/(6m_1^3)$", "(b)"),
              ("D_left", "SE_D_left", 0.0, r"$D$ (left endpoint)", "(c)")]

    for ax, (name, sname, target, lab, tag) in zip(axes, panels):
        y, s = loc[name].to_numpy(), loc[sname].to_numpy()
        r = plr.loc[name]
        fit = target + np.sign(y[-1] - target) * np.exp(r.logA) * Lf ** (-r.beta)
        ax.plot(Lf, fit, "-", color=C["red"], lw=1.0, zorder=2)
        ax.errorbar(L, y, yerr=s, fmt="o", ms=3, capsize=1.5, lw=0.8,
                    color=C["blue"], zorder=3)
        ax.set_xlabel(r"$\log X_{\mathrm{eff}}$")
        ax.set_ylabel(lab)
        _panel(ax, tag, rf"$\beta={r.beta:.3f}$")
        _pad_ylim(ax, np.r_[y - s, y + s, fit])

    handles = [Line2D([], [], color=C["blue"], marker="o", ms=3, ls="none",
                      label=r"census ($\pm 1$ bootstrap SE)"),
               Line2D([], [], color=C["red"], lw=1.2,
                      label=r"fit $\propto(\log X)^{-\beta}$")]
    _below(fig, handles, [h.get_label() for h in handles], ncol=2)
    _save(fig, outdir, "fig3_local")


# ------------------------------------------------------------ Figure 4
def fig_model(mvd: pd.DataFrame, outdir):
    """R(g;N) against the parameter-free density model, and the residual.

    Colour carries the scale and line style carries the source, so the legend
    factorises into three scales plus three sources instead of nine entries.
    """
    fig, axes = plt.subplots(1, 2, figsize=(FULL, 2.6), layout="constrained")
    cols = [C["blue"], C["orange"], C["green"]]
    ks = sorted(mvd.k.unique())
    lab_of = {}

    for i, k in enumerate(ks):
        grp = mvd[mvd.k == k].sort_values("g")
        lab_of[k] = rf"$X=10^{{{math.log10(grp.X.iloc[0]):.1f}}}$"
        axes[0].plot(grp.g, grp.R_log, "-", lw=1.0, color=cols[i], zorder=2)
        axes[0].plot(grp.g, grp.R_pi, "--", lw=0.9, color=cols[i], zorder=2)
        axes[0].plot(grp.g, grp.R, "o", ms=2.6, color=cols[i], zorder=3)
        axes[1].plot(grp.g, grp.R - grp.R_log, "o-", ms=2.4, lw=0.8,
                     color=cols[i], zorder=3)

    axes[0].axhline(1, color=C["black"], lw=0.7, ls=":", zorder=1)
    axes[0].set_xlabel("gap $g$")
    axes[0].set_ylabel("$R(g;N)$")
    _panel(axes[0], "(a)", "census against model")
    _pad_ylim(axes[0], np.r_[mvd.R.to_numpy(), mvd.R_log.to_numpy(),
                             mvd.R_pi.to_numpy()])

    axes[1].axhline(0, color=C["black"], lw=0.7, ls=":", zorder=1)
    axes[1].set_xlabel("gap $g$")
    axes[1].set_ylabel(r"$R_{\mathrm{census}} - R_{\mathrm{model}}$")
    _panel(axes[1], "(b)", r"residual, $1/\log t$ envelope")
    _pad_ylim(axes[1], (mvd.R - mvd.R_log).to_numpy())

    scale_h = [Line2D([], [], color=cols[i], lw=1.4, label=lab_of[k])
               for i, k in enumerate(ks)]
    src_h = [Line2D([], [], color=C["grey"], marker="o", ms=3, ls="none",
                    label="census"),
             Line2D([], [], color=C["grey"], lw=1.4, ls="-",
                    label=r"model, $1/\log t$"),
             Line2D([], [], color=C["grey"], lw=1.2, ls="--",
                    label=r"model, $\pi(t)/t$")]
    _below(fig, scale_h + src_h, [h.get_label() for h in scale_h + src_h], ncol=6)
    _save(fig, outdir, "fig4_model")


# ------------------------------------------------------------ Figure 5
def fig_collapse(spat: pd.DataFrame, outdir, min_count=2000, xmax=2.6):
    """The scaling collapse. Scale is a sequential variable, so it is encoded
    by a colourbar rather than by a legend naming three of the twenty curves."""
    ks = sorted(spat.k.unique())
    lo = math.log10(spat[spat.k == ks[0]].X.iloc[0])
    hi = math.log10(spat[spat.k == ks[-1]].X.iloc[0])
    norm = matplotlib.colors.Normalize(vmin=lo, vmax=hi)
    cmap = plt.get_cmap(SEQ_CMAP)

    fig, axes = plt.subplots(1, 3, figsize=(FULL, 2.5), layout="constrained")
    ymain, ycoll = [], []

    for k in ks:
        grp = spat[(spat.k == k) & (spat["count"] >= min_count)
                   & (spat.x <= xmax)].sort_values("g")
        if grp.empty:
            continue
        col = cmap(norm(math.log10(grp.X.iloc[0])))
        L = grp.logX.iloc[0]
        axes[0].plot(grp.g, grp.R, "-", lw=0.7, color=col)
        axes[1].plot(grp.x, grp.R, "-", lw=0.7, color=col)
        axes[2].plot(grp.x, L * (grp.R - 1), "-", lw=0.7, color=col)
        ymain.append(grp.R.to_numpy())
        ycoll.append((L * (grp.R - 1)).to_numpy())

    for ax in axes[:2]:
        ax.axhline(1, color=C["black"], lw=0.7, ls=":", zorder=0)
        ax.set_ylabel("$R(g;N)$")
        _pad_ylim(ax, np.concatenate(ymain))
    axes[0].set_xlabel("gap $g$")
    axes[1].axvline(1, color=C["black"], lw=0.7, ls="--", zorder=0)
    axes[1].set_xlabel(r"$x = g/\log X$")
    axes[2].axhline(0, color=C["black"], lw=0.7, ls=":", zorder=0)
    axes[2].axvline(1, color=C["black"], lw=0.7, ls="--", zorder=0)
    axes[2].set_xlabel(r"$x = g/\log X$")
    axes[2].set_ylabel(r"$\log X\,(R(g;N)-1)$")

    _panel(axes[0], "(a)", "unscaled")
    _panel(axes[1], "(b)", "abscissa rescaled")
    _panel(axes[2], "(c)", "both rescaled")

    xx = np.linspace(0.08, xmax, 60)
    Lmax = float(spat.logX.max())
    m_log = [Lmax * (R_model(x * Lmax, Lmax, "log") - 1) for x in xx]
    m_pi = [Lmax * (R_model(x * Lmax, Lmax, "pi") - 1) for x in xx]
    first = (xx - 1) / 2
    axes[2].plot(xx, m_log, "-", color=C["red"], lw=1.3, zorder=5)
    axes[2].plot(xx, m_pi, "-.", color=C["orange"], lw=1.1, zorder=5)
    axes[2].plot(xx, first, "--", color=C["grey"], lw=1.0, zorder=5)
    _pad_ylim(axes[2], np.concatenate([np.concatenate(ycoll), m_log, m_pi, first]))

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=axes.tolist(), location="right", fraction=0.021,
                      pad=0.012, aspect=28)
    cb.set_label(r"$\log_{10} X$", fontsize=7.5)
    cb.ax.tick_params(labelsize=6.5, width=0.5)
    cb.outline.set_linewidth(0.5)

    handles = [Line2D([], [], color=C["red"], lw=1.4,
                      label=rf"density model, $\log X={Lmax:.0f}$"),
               Line2D([], [], color=C["orange"], lw=1.2, ls="-.",
                      label=r"model, $\pi(t)/t$ envelope"),
               Line2D([], [], color=C["grey"], lw=1.1, ls="--",
                      label=r"first order, $(x-1)/2$")]
    _below(fig, handles, [h.get_label() for h in handles], ncol=3)
    _save(fig, outdir, "fig5_collapse")


# ------------------------------------------------------------ Figure 6
def fig_crossover(cross: pd.DataFrame, asym: pd.DataFrame, outdir):
    """The crossover under several locators, and the deficit against the model.

    Panel (b) is drawn over the range the census covers. The model's asymptote
    is shown as a horizontal reference instead of extending the axis to log X =
    50, which previously left most of the panel empty.
    """
    fig, axes = plt.subplots(1, 2, figsize=(FULL, 2.7), layout="constrained")
    L = cross.logX.to_numpy()

    axes[0].plot(L, L, ":", color=C["black"], lw=0.9, zorder=1)
    axes[0].plot(L, L - 1.5, "--", color=C["red"], lw=0.9, zorder=1)
    axes[0].plot(L, cross.gstar_x515, "s", ms=2.6, color=C["sky"], zorder=2)
    axes[0].plot(L, cross.gstar_b02, "^", ms=2.6, color=C["orange"], zorder=2)
    axes[0].plot(L, cross.gstar_b05, "v", ms=2.6, color=C["green"], zorder=2)
    axes[0].errorbar(L, cross.gstar_wls, yerr=cross.SE_gstar, fmt="o", ms=3.2,
                     lw=0.8, capsize=1.5, color=C["blue"], zorder=3)
    axes[0].set_xlabel(r"$\log X$")
    axes[0].set_ylabel(r"$g^*(X)$")
    _panel(axes[0], "(a)", "crossover under six locators")
    _pad_ylim(axes[0], np.r_[L, L - 1.5, cross.gstar_wls - cross.SE_gstar,
                             cross.gstar_wls + cross.SE_gstar,
                             cross.gstar_x515.to_numpy(),
                             cross.gstar_b02.to_numpy(),
                             cross.gstar_b05.to_numpy()])

    locs = ["gstar_x713", "gstar_x614", "gstar_x515",
            "gstar_b02", "gstar_b03", "gstar_b05"]
    allc = cross[locs]
    dfc = L - cross.gstar_wls.to_numpy()
    se = cross.SE_gstar.to_numpy()
    xlo, xhi = L.min() - 0.5, L.max() + 0.5

    # The model is tabulated on a sparse, very wide grid (log X up to 230), so
    # restricting it to the census range leaves too few points to draw. Its
    # deficit behaves as c0 + c1/log X, so interpolate in 1/log X, where it is
    # close to linear, and evaluate on a dense grid over the census range.
    a = asym.sort_values("logX")
    grid = np.linspace(xlo, xhi, 120)
    inv = 1.0 / a.logX.to_numpy()
    order = np.argsort(inv)
    model_curve = np.interp(1.0 / grid, inv[order], a.deficit.to_numpy()[order])

    axes[1].fill_between(L, L - allc.max(axis=1), L - allc.min(axis=1),
                         color=C["sky"], alpha=0.30, lw=0, zorder=1)
    axes[1].axhline(1.5, color=C["red"], lw=0.9, ls="--", zorder=2)
    axes[1].plot(grid, model_curve, "-", color=C["red"], lw=1.2, zorder=3)
    axes[1].errorbar(L, dfc, yerr=se, fmt="o", ms=3, lw=0.8, capsize=1.5,
                     color=C["blue"], zorder=4)
    axes[1].set_xlim(xlo, xhi)
    axes[1].set_xlabel(r"$\log X$")
    axes[1].set_ylabel(r"$\log X - g^*$")
    _panel(axes[1], "(b)", "deficit against the model")
    _pad_ylim(axes[1], np.r_[dfc - se, dfc + se, 1.5, model_curve,
                             (L - allc.max(axis=1)).to_numpy(),
                             (L - allc.min(axis=1)).to_numpy()])

    handles = [
        Line2D([], [], color=C["blue"], marker="o", ms=3.2, ls="none",
               label=r"WLS, $x\in[0.6,1.4]$ (primary, $\pm 1$ SE)"),
        Line2D([], [], color=C["sky"], marker="s", ms=2.8, ls="none",
               label=r"window $x\in[0.5,1.5]$"),
        Line2D([], [], color=C["orange"], marker="^", ms=2.8, ls="none",
               label=r"band $|R-1|<0.02$"),
        Line2D([], [], color=C["green"], marker="v", ms=2.8, ls="none",
               label=r"band $|R-1|<0.05$"),
        Line2D([], [], color=C["black"], lw=1.0, ls=":", label=r"$g^*=\log X$"),
        Line2D([], [], color=C["red"], lw=1.2, label="density model"),
        Line2D([], [], color=C["red"], lw=1.0, ls="--",
               label=r"model asymptote $3/2$"),
        Patch(facecolor=C["sky"], alpha=0.30, lw=0,
              label="range over six locators"),
    ]
    _below(fig, handles, [h.get_label() for h in handles], ncol=4)
    _save(fig, outdir, "fig6_crossover")
