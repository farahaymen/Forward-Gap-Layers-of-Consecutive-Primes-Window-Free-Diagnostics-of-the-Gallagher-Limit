#!/usr/bin/env python3
"""
run_all.py -- regenerate every number, table and figure in the paper.

    python3 run_all.py ../data/fgp_layers.csv --out ../paper --verify-limit 2100000000

Outputs (under --out):
    tables/*.tex      LaTeX tabulars, included by main.tex
    tables/macros.tex numbers quoted in the running text
    figures/*.pdf     figures 1-6 (and PNG previews)
    csv/*.csv         every intermediate table
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fgp import core, frequency, local, spatial, figures, tables, verify  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("census")
    ap.add_argument("--out", default="../paper")
    ap.add_argument("--verify-limit", type=int, default=2_100_000_000,
                    help="upper bound X for the independent NumPy census (0 = skip)")
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--sieve-seconds", type=float, default=14.3,
                    help="wall time of fgp_sieve.c on the reference machine (quoted only)")
    a = ap.parse_args()

    tdir, fdir, cdir = (os.path.join(a.out, d) for d in ("tables", "figures", "csv"))
    for d in (tdir, fdir, cdir):
        os.makedirs(d, exist_ok=True)

    t0 = time.time()
    scales = core.load_census(a.census)
    print(f"loaded {len(scales)} scales, X from {scales[0].X:,} to {scales[-1].X:,}")

    # ---- Section 3: fitted rate and specification grid
    grid = frequency.specification_grid(scales)
    prim = frequency.primary_table(scales, grid)
    summ = frequency.grid_summary(grid)
    wolf = frequency.wolf_comparison(grid, prim)
    grid.to_csv(f"{cdir}/specification_grid.csv", index=False)
    prim.to_csv(f"{cdir}/census_primary.csv", index=False)
    summ.to_csv(f"{cdir}/specification_summary.csv", index=False)
    print(f"Wolf s={wolf['published']} at X~{wolf['X']:,}: ours {wolf['s_primary']:.4f}, "
          f"envelope [{wolf['s_min']:.4f},{wolf['s_max']:.4f}] -> {'INSIDE' if wolf['inside'] else 'outside'}; "
          f"closest spec g>={wolf['best_gmin']}, c={wolf['best_c']}, {wolf['best_est']} gives {wolf['best_s']:.4f}")
    print(f"C_fit primary: {prim.C_fit.iloc[0]:.4f} -> {prim.C_fit.iloc[-1]:.4f}; "
          f"{int(summ.monotone_dec.sum())} monotone decreasing, {int(summ.monotone_inc.sum())} increasing")

    # ---- Section 4: local diagnostics
    loc = local.local_table(scales, B=a.bootstrap)
    st = local.step_tests(loc)
    pl = local.power_law_fits(loc)
    loc.to_csv(f"{cdir}/local_diagnostics.csv", index=False)
    st.drop(columns="z_steps").to_csv(f"{cdir}/step_tests.csv", index=False)
    pl.to_csv(f"{cdir}/power_law_fits.csv", index=False)
    print(f"E2: {loc.E2.iloc[0]:.5f} -> {loc.E2.iloc[-1]:.5f}; steps toward limit: "
          + ", ".join(f"{r.statistic} {r.n_toward}/{r.n_steps}" for _, r in st.iterrows()))
    print("power-law exponents: " + ", ".join(f"{r.statistic} {r.beta:.3f}+-{r.SE_beta:.3f}" for _, r in pl.iterrows()))

    # ---- Section 5: spatial statistic, crossover, model
    spat = spatial.spatial_table(scales)
    cross = spatial.crossover_table(spat)
    asym = spatial.model_asymptotics("log")
    asym_pi = spatial.model_asymptotics("pi")
    ks3 = [scales[0].k, scales[9].k, scales[-1].k]
    mvd = spatial.model_vs_data(spat, ks3)
    amp = spatial.amplitude_ratio(spat)
    amp.to_csv(f'{cdir}/amplitude_ratio.csv', index=False)
    print(f"amplitude ratio (census/model slope): {amp.ratio_log.mean():.3f} +- {amp.ratio_log.std(ddof=1):.3f} "
          f"over {len(amp)} scales; pi-variant {amp.ratio_pi.mean():.3f}")
    spat.to_csv(f"{cdir}/spatial.csv", index=False)
    cross.to_csv(f"{cdir}/crossover.csv", index=False)
    asym.to_csv(f"{cdir}/model_asymptotics_log.csv", index=False)
    asym_pi.to_csv(f"{cdir}/model_asymptotics_pi.csv", index=False)
    mvd.to_csv(f"{cdir}/model_vs_data.csv", index=False)
    print(f"g*: {cross.gstar.iloc[0]:.2f} -> {cross.gstar.iloc[-1]:.2f}; g*/logX mean {cross.ratio.mean():.3f} sd {cross.ratio.std(ddof=1):.3f}")
    print(f"model deficit at log X = 23: {asym[asym.log10X==10].deficit.iloc[0]:.4f}; Richardson -> {asym.attrs['richardson_c0']:.4f}")
    print(f"pi-variant model deficit at log X = 23: {asym_pi[asym_pi.log10X==10].deficit.iloc[0]:.4f}")
    for k, grp in mvd.groupby("k"):
        print(f"  k={k}: max |R - R_log| for g<=30: {grp[grp.g<=30].eval('abs(R-R_log)').max():.4f}, "
              f"pi-variant: {grp[grp.g<=30].eval('abs(R-R_pi)').max():.4f}")

    # ---- Section 6: independent verification
    if a.verify_limit > 0:
        ver = verify.independent_census(scales, a.verify_limit)
        ver.to_csv(f"{cdir}/verification.csv", index=False)
        print(f"independent census to {ver.X.max():,}: {len(ver)} layer counts compared, "
              f"{int((~ver.match).sum())} mismatches, {ver.attrs['seconds']:.0f} s")
    else:
        import pandas as pd
        ver = pd.DataFrame(dict(k=[], X=[], g=[], count_C=[], count_py=[], match=[]))

    # ---- tables, macros, figures
    tables.table_census(prim, tdir)
    tables.table_specgrid(summ, tdir)
    tables.table_local(loc, tdir)
    tables.table_uncertainty(loc, tdir)
    tables.table_steps(st, pl, tdir)
    tables.table_crossover(cross, tdir)
    tables.table_verify(ver, tdir)
    tables.write_macros(prim, loc, st, pl, cross, asym, summ, ver, a.sieve_seconds, tdir, wolf, amp)

    figures.fig_cfit(prim, fdir)
    figures.fig_specgrid(grid, fdir)
    figures.fig_local(loc, pl, fdir)
    figures.fig_model(mvd, fdir)
    figures.fig_collapse(spat, fdir)
    figures.fig_crossover(cross, asym, fdir)
    print(f"done in {time.time() - t0:.0f} s; outputs under {a.out}")


if __name__ == "__main__":
    main()
