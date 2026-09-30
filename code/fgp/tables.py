"""
fgp.tables -- write the paper's tables as LaTeX fragments (booktabs).

Each function writes a file containing only the tabular environment, so the
manuscript wraps it in \\begin{table} with its own caption.
"""

import os

import numpy as np
import pandas as pd


def _w(outdir, name, body):
    with open(os.path.join(outdir, name + ".tex"), "w") as fh:
        fh.write(body)


def _n(x):
    return f"{int(x):,}".replace(",", "\\,")


def table_census(prim: pd.DataFrame, outdir):
    lines = [r"\begin{tabular}{rrrcccccc}", r"\toprule",
             r"$N$ & $X=p_N$ & $\log X$ & $C_W(X)$ & $C_{\mathrm{fit}}(X)$ & SE & $\Delta(X)$ & envelope & pts \\",
             r"\midrule"]
    for _, r in prim.iterrows():
        lines.append(f"{_n(r.N)} & {_n(r.X)} & {r.logX:.3f} & {r.C_W:.4f} & {r.C_fit:.4f} & "
                     f"{r.SE:.5f} & {r.delta:.4f} & [{r.env_min:.3f}, {r.env_max:.3f}] & {int(r.n_points)} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    _w(outdir, "tab_census", "\n".join(lines))


def table_specgrid(summ: pd.DataFrame, outdir):
    """The specification grid, split into two side-by-side blocks by g_min so
    that it occupies half a page rather than a full one."""
    est_lab = {"POIS": "Poisson", "WLS": "WLS", "OLS": "OLS"}
    blocks = []
    for gmin in sorted(summ.g_min.unique()):
        b = summ[summ.g_min == gmin].sort_values(["c", "estimator"]).reset_index(drop=True)
        blocks.append(b)
    n = max(len(b) for b in blocks)
    head = (r"window & est. & first & last & net & mon.")
    lines = [r"\begin{tabular}{@{}cl rrr c@{\qquad} cl rrr c@{}}", r"\toprule",
             r"\multicolumn{6}{c}{$g \geq " + str(int(blocks[0].g_min.iloc[0])) + r"$} & "
             r"\multicolumn{6}{c}{$g \geq " + str(int(blocks[1].g_min.iloc[0])) + r"$} \\",
             r"\cmidrule(r){1-6}\cmidrule(l){7-12}",
             head + " & " + head + r" \\", r"\midrule"]
    for i in range(n):
        cells = []
        for b in blocks:
            if i < len(b):
                r = b.iloc[i]
                mono = "dec." if r.monotone_dec else ("inc." if r.monotone_inc else "--")
                cells.append(f"${int(r.c)}\\log X$ & {est_lab[r.estimator]} & "
                             f"{r.C_first:.4f} & {r.C_last:.4f} & {r.net:+.4f} & {mono}")
            else:
                cells.append(" & & & & & ")
        lines.append(" & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    _w(outdir, "tab_specgrid", "\n".join(lines))


def table_local(loc: pd.DataFrame, outdir):
    lines = [r"\begin{tabular}{rrr c c c c c c}", r"\toprule",
             r"$k$ & $X_{k-1}$ & $X_k$ & gaps & $m_1$ & $\log X_{\mathrm{eff}}/m_1$ & $E_2$ & $E_3$ & $D$ \\",
             r"\midrule"]
    for _, r in loc.iterrows():
        lines.append(f"{int(r.k)} & {_n(r.X_lo)} & {_n(r.X_hi)} & {_n(r.n_gaps)} & {r.m1:.4f} & "
                     f"{r.ratio:.5f} & {r.E2:.5f} & {r.E3:.5f} & {r.D_left:.5f} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    _w(outdir, "tab_local", "\n".join(lines))


def table_uncertainty(loc: pd.DataFrame, outdir, ks=(1, 5, 9, 13, 16, 19)):
    lines = [r"\begin{tabular}{r c c c c c}", r"\toprule",
             r"$k$ & $\log X_{\mathrm{eff}}$ & $E_2$ (SE) & $E_3$ (SE) & $D_{\mathrm{left}}$ (SE) & $D_{\mathrm{mid}}$ (SE) \\",
             r"\midrule"]
    for _, r in loc[loc.k.isin(ks)].iterrows():
        lines.append(f"{int(r.k)} & {r.logX_eff:.3f} & {r.E2:.5f} ({r.SE_E2:.5f}) & {r.E3:.5f} ({r.SE_E3:.5f}) & "
                     f"{r.D_left:.5f} ({r.SE_D_left:.5f}) & {r.D_mid:.5f} ({r.SE_D_mid:.5f}) \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    _w(outdir, "tab_uncertainty", "\n".join(lines))


def table_steps(st: pd.DataFrame, pl: pd.DataFrame, outdir):
    plr = pl.set_index("statistic")
    names = {"E2": "$E_2$", "E3": "$E_3$", "D_left": "$D_{\\mathrm{left}}$",
             "D_mid": "$D_{\\mathrm{mid}}$", "W2": "$W^2$"}
    lines = [r"\begin{tabular}{l c c c c c c}", r"\toprule",
             r"statistic & steps toward limit & median $z$ & min $z$ & end-to-end change & end-to-end $z$ & exponent $\beta$ \\",
             r"\midrule"]
    for _, r in st.iterrows():
        beta = f"{plr.loc[r.statistic].beta:.3f} $\\pm$ {plr.loc[r.statistic].SE_beta:.3f}" \
            if r.statistic in plr.index else "--"
        lines.append(f"{names[r.statistic]} & {int(r.n_toward)}/{int(r.n_steps)} & {r.median_z:.1f} & "
                     f"{r.min_z:.2f} & {r.end_change:+.4f} $\\pm$ {r.end_SE:.4f} & {r.end_z:.0f} & {beta} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    _w(outdir, "tab_steps", "\n".join(lines))


def table_crossover(cross: pd.DataFrame, outdir):
    lines = [r"\begin{tabular}{rr ccc ccc c c}", r"\toprule",
             r"$X=p_N$ & $\log X$ & \multicolumn{3}{c}{window in $x$} & \multicolumn{3}{c}{band in $R$} & \multicolumn{2}{c}{weighted fit} \\",
             r"\cmidrule(lr){3-5}\cmidrule(lr){6-8}",
             r" & & $[0.7,1.3]$ & $[0.6,1.4]$ & $[0.5,1.5]$ & 0.02 & 0.03 & 0.05 & $g^*$ (WLS) & $\log X-g^*$ \\",
             r"\midrule"]
    for _, r in cross.iterrows():
        lines.append(f"{_n(r.X)} & {r.logX:.3f} & {r.gstar_x713:.2f} & \\textbf{{{r.gstar_x614:.2f}}} & {r.gstar_x515:.2f} & "
                     f"{r.gstar_b02:.2f} & {r.gstar_b03:.2f} & {r.gstar_b05:.2f} & "
                     f"{r.gstar_wls:.2f} $\\pm$ {r.SE_gstar:.2f} & {r.logX - r.gstar_wls:.2f} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    _w(outdir, "tab_crossover", "\n".join(lines))


def table_verify(ver: pd.DataFrame, outdir):
    g = ver.groupby("k").agg(X=("X", "first"), rows=("g", "size"), mismatches=("match", lambda m: int((~m).sum())))
    lines = [r"\begin{tabular}{r r r r}", r"\toprule",
             r"checkpoint $k$ & $X$ & layers compared & mismatches \\", r"\midrule"]
    for k, r in g.iterrows():
        lines.append(f"{int(k)} & {_n(r.X)} & {int(r.rows)} & {int(r.mismatches)} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    _w(outdir, "tab_verify", "\n".join(lines))


def write_macros(prim, loc, st, pl, cross, asym, summ, ver, seconds, outdir, wolf=None, amp=None):
    """Numbers quoted in the running text, as LaTeX macros."""
    plr = pl.set_index("statistic")
    str_ = st.set_index("statistic")
    n_dec = int(summ.monotone_dec.sum())
    n_inc = int(summ.monotone_inc.sum())
    k9 = prim[prim.N == 9_000_000].iloc[0]
    m = {
        "nScales": len(prim),
        "Nfirst": _n(prim.N.iloc[0]), "Nlast": _n(prim.N.iloc[-1]),
        "Xfirst": _n(prim.X.iloc[0]), "Xlast": _n(prim.X.iloc[-1]),
        "CfitFirst": f"{prim.C_fit.iloc[0]:.4f}", "CfitLast": f"{prim.C_fit.iloc[-1]:.4f}",
        "CWFirst": f"{prim.C_W.iloc[0]:.4f}", "CWLast": f"{prim.C_W.iloc[-1]:.4f}",
        "deltaFirst": f"{prim.delta.iloc[0]:.4f}", "deltaLast": f"{prim.delta.iloc[-1]:.4f}",
        "spreadMin": f"{prim.env_spread.min():.3f}", "spreadMax": f"{prim.env_spread.max():.3f}",
        "netChangePrimary": f"{prim.C_fit.iloc[0] - prim.C_fit.iloc[-1]:.3f}",
        "nMonoDec": n_dec, "nMonoInc": n_inc,
        "envNineMin": f"{k9.env_min:.4f}", "envNineMax": f"{k9.env_max:.4f}",
        "narrowLast": f"{summ.C_last.min():.4f}",
        "envMaxLast": f"{summ.C_last.max():.4f}",
        "EtwoFirst": f"{loc.E2.iloc[0]:.5f}", "EtwoLast": f"{loc.E2.iloc[-1]:.5f}",
        "EthreeFirst": f"{loc.E3.iloc[0]:.5f}", "EthreeLast": f"{loc.E3.iloc[-1]:.5f}",
        "DFirst": f"{loc.D_left.iloc[0]:.5f}", "DLast": f"{loc.D_left.iloc[-1]:.5f}",
        "DmidFirst": f"{loc.D_mid.iloc[0]:.5f}", "DmidLast": f"{loc.D_mid.iloc[-1]:.5f}",
        "ratioMin": f"{loc.ratio.min():.5f}", "ratioMax": f"{loc.ratio.max():.5f}",
        "betaEtwo": f"{plr.loc['E2'].beta:.3f}", "betaEtwoSE": f"{plr.loc['E2'].SE_beta:.3f}",
        "betaEthree": f"{plr.loc['E3'].beta:.3f}", "betaEthreeSE": f"{plr.loc['E3'].SE_beta:.3f}",
        "betaD": f"{plr.loc['D_left'].beta:.3f}", "betaDSE": f"{plr.loc['D_left'].SE_beta:.3f}",
        "EtwoAtFifty": f"{plr.loc['E2'].at_50:.3f}", "EtwoAtTwoHundred": f"{plr.loc['E2'].at_200:.3f}",
        "zEtwoEnd": f"{str_.loc['E2'].end_z:.0f}", "zEthreeEnd": f"{str_.loc['E3'].end_z:.0f}",
        "zDEnd": f"{str_.loc['D_left'].end_z:.0f}",
        "medzEtwo": f"{str_.loc['E2'].median_z:.1f}", "medzEthree": f"{str_.loc['E3'].median_z:.1f}",
        "medzD": f"{str_.loc['D_left'].median_z:.1f}",
        "towardEtwo": int(str_.loc['E2'].n_toward), "towardEthree": int(str_.loc['E3'].n_toward),
        "towardD": int(str_.loc['D_left'].n_toward), "nSteps": int(str_.loc['E2'].n_steps),
        "zFirstStepEtwo": f"{str_.loc['E2'].z_steps[0]:.2f}",
        "zFirstStepEthree": f"{str_.loc['E3'].z_steps[0]:.2f}",
        "zFirstStepD": f"{str_.loc['D_left'].z_steps[0]:.2f}",
        "gstarFirst": f"{cross.gstar.iloc[0]:.2f}", "gstarLast": f"{cross.gstar.iloc[-1]:.2f}",
        "ratioMean": f"{cross.ratio.mean():.3f}", "ratioSD": f"{cross.ratio.std(ddof=1):.3f}",
        "deficitFirst": f"{cross.deficit.iloc[0]:.2f}", "deficitLast": f"{cross.deficit.iloc[-1]:.2f}",
        "bandSpreadMax": f"{(cross[['gstar_x713','gstar_x614','gstar_x515','gstar_b02','gstar_b03','gstar_b05']].max(axis=1) - cross[['gstar_x713','gstar_x614','gstar_x515','gstar_b02','gstar_b03','gstar_b05']].min(axis=1)).max():.2f}",
        "ratioAllMin": f"{cross[['gstar_x713','gstar_x614','gstar_x515','gstar_b02','gstar_b03','gstar_b05']].div(cross.logX, axis=0).min().min():.3f}",
        "ratioAllMax": f"{cross[['gstar_x713','gstar_x614','gstar_x515','gstar_b02','gstar_b03','gstar_b05']].div(cross.logX, axis=0).max().max():.3f}",
        "deficitMin": f"{cross.deficit.min():.2f}", "deficitMax": f"{cross.deficit.max():.2f}",
        "deficitMean": f"{cross.deficit.mean():.2f}", "deficitSD": f"{cross.deficit.std(ddof=1):.2f}",
        "gstarLastSE": f"{cross.SE_gstar.iloc[-1]:.2f}",
        "deficitLastSE": f"{cross.SE_gstar.iloc[-1]:.2f}",
        "deficitLastVal": f"{(cross.logX - cross.gstar_wls).iloc[-1]:.2f}",
        "ratioLastVal": f"{(cross.gstar_wls / cross.logX).iloc[-1]:.4f}",
        "ratioLastSE": f"{(cross.SE_gstar / cross.logX).iloc[-1]:.4f}",
        "sigmaUnity": f"{((cross.logX - cross.gstar_wls) / cross.SE_gstar).iloc[-1]:.1f}",
        "sigmaModel": f"{(abs((cross.logX - cross.gstar_wls) - 1.5) / cross.SE_gstar).iloc[-1]:.1f}",
        "modelDeficitTen": f"{asym[asym.log10X == 10].deficit.iloc[0]:.4f}",
        "modelDeficitHundred": f"{asym[asym.log10X == 100].deficit.iloc[0]:.4f}",
        "richardson": f"{asym.attrs['richardson_c0']:.4f}",
        "verifyLayers": _n(len(ver)), "verifyMismatch": int((~ver.match.astype(bool)).sum()) if len(ver) else 0,
        "verifyXmax": _n(ver.X.max()) if len(ver) else "0", "verifySeconds": f"{ver.attrs.get('seconds', 0):.0f}",
        "sieveSeconds": f"{seconds:.0f}",
        "nSpecs": len(summ),
        "sFirst": f"{prim.s_wolf.iloc[0]:.4f}", "sLast": f"{prim.s_wolf.iloc[-1]:.4f}",
    }
    if amp is not None:
        m.update({
            "ampRatio": f"{amp.ratio_log.mean():.2f}",
            "ampRatioSD": f"{amp.ratio_log.std(ddof=1):.2f}",
            "ampRatioMin": f"{amp.ratio_log.min():.2f}",
            "ampRatioMax": f"{amp.ratio_log.max():.2f}",
            "ampRatioPi": f"{amp.ratio_pi.mean():.2f}",
            "ampRatioPiSD": f"{amp.ratio_pi.std(ddof=1):.2f}",
            "ampSlopeData": f"{amp.slope_data.iloc[-1]:.2f}",
            "ampSlopeModel": f"{amp.slope_log.iloc[-1]:.2f}",
        })
    if wolf is not None:
        est_name = {"POIS": "Poisson GLM", "WLS": "weighted least squares",
                    "OLS": "unweighted least squares"}[wolf["best_est"]]
        m.update({
            "wolfX": _n(wolf["X"]),
            "wolfPublished": f"{wolf['published']:.3f}",
            "wolfOurs": f"{wolf['s_primary']:.4f}",
            "wolfEnvMin": f"{wolf['s_min']:.4f}", "wolfEnvMax": f"{wolf['s_max']:.4f}",
            "wolfInside": "inside" if wolf["inside"] else "outside",
            "wolfBestC": wolf["best_c"], "wolfBestGmin": wolf["best_gmin"],
            "wolfBestEst": est_name, "wolfBestS": f"{wolf['best_s']:.4f}",
        })
    with open(os.path.join(outdir, "macros.tex"), "w") as fh:
        for k, v in m.items():
            fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
