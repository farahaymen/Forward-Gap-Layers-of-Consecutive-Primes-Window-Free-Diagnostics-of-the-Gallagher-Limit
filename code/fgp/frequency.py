"""
fgp.frequency -- the fitted rescaled decay rate C_fit(X) and its specification grid.

Model:  E|L(g;N)| = exp(a + log H(g) - lambda * g), fitted over even g in a window.
C_fit(X) = lambda * log X.  The Wolf benchmark is C_W(X) = pi(X) log X / X,
which equals the method-of-moments estimator of the same rate (mean gap = X/pi(X)).

Three estimators:
    POIS  Poisson log-linear regression by IRLS, log H(g) as an offset
    WLS   weighted least squares of log(f/H) on g, weights = counts
    OLS   unweighted least squares of log(f/H) on g
Windows: g_min in {2, 4}, g_max = c * log X for c in {1, 2, 3, 4}.
"""

import math
from itertools import product

import numpy as np
import pandas as pd

from .core import H, Scale

ESTIMATORS = ("POIS", "WLS", "OLS")
GMINS = (2, 4)
# Windows g <= c log X. c runs to 6 because Wolf's published fit, described as
# taken over the 'linear portions' of a log plot, corresponds to an unweighted
# fit with c near 5; a grid stopping at c = 4 would not cover the literature.
CS = (1, 2, 3, 4, 5, 6)
PRIMARY = (4, 3, "POIS")


def fit_poisson(gs, cs):
    """IRLS for log mu = a + log H(g) - lam*g. Returns (lam, SE(lam))."""
    x = np.asarray(gs, float)
    y = np.asarray(cs, float)
    off = np.array([math.log(H(int(g))) for g in gs])
    Xd = np.column_stack([np.ones_like(x), x])
    b = np.array([math.log(max(y.mean(), 1e-9)), -0.05])
    for _ in range(200):
        eta = Xd @ b + off
        mu = np.exp(np.clip(eta, -700, 700))
        z = eta - off + (y - mu) / np.maximum(mu, 1e-12)
        A = Xd.T @ (mu[:, None] * Xd)
        bn = np.linalg.solve(A, Xd.T @ (mu * z))
        done = np.max(np.abs(bn - b)) < 1e-12
        b = bn
        if done:
            break
    mu = np.exp(np.clip(Xd @ b + off, -700, 700))
    cov = np.linalg.inv(Xd.T @ (mu[:, None] * Xd))
    return -b[1], math.sqrt(cov[1, 1])


def fit_ls(gs, cs, tot, weighted):
    x = np.asarray(gs, float)
    y = np.array([math.log(c / tot / H(int(g))) for g, c in zip(gs, cs)])
    w = np.asarray(cs, float) if weighted else np.ones_like(x)
    Xd = np.column_stack([np.ones_like(x), x])
    A = Xd.T @ (w[:, None] * Xd)
    b = np.linalg.solve(A, Xd.T @ (w * y))
    r = y - Xd @ b
    s2 = float((w * r * r).sum() / (len(x) - 2))
    cov = s2 * np.linalg.inv(A)
    return -b[1], math.sqrt(cov[1, 1])


def specification_grid(scales: list[Scale]) -> pd.DataFrame:
    rows = []
    for s in scales:
        logX, tot = s.logX, s.N - 1
        for gmin, c, est in product(GMINS, CS, ESTIMATORS):
            gmax = int(c * logX) // 2 * 2
            mask = (s.g >= gmin) & (s.g <= gmax) & (s.count > 0)
            gs, cs = s.g[mask], s.count[mask]
            if len(gs) < 4:
                continue
            if est == "POIS":
                lam, se = fit_poisson(gs, cs)
            else:
                lam, se = fit_ls(gs, cs, tot, est == "WLS")
            rows.append(dict(k=s.k, N=s.N, X=s.X, logX=logX, g_min=gmin, c=c,
                             estimator=est, g_max=gmax, n_points=int(mask.sum()),
                             C_fit=lam * logX, SE=se * logX))
    return pd.DataFrame(rows)


def primary_table(scales: list[Scale], grid: pd.DataFrame) -> pd.DataFrame:
    gmin, c, est = PRIMARY
    prim = grid[(grid.g_min == gmin) & (grid.c == c) & (grid.estimator == est)]
    prim = prim.set_index("k")
    rows = []
    for s in scales:
        r = prim.loc[s.k]
        env = grid[grid.k == s.k].C_fit
        CW = s.N * s.logX / s.X
        # Wolf fits a(x)exp(-s u) in the variable u = g*pi(X)/X, so his slope
        # satisfies lambda = s*pi(X)/X and hence s = C_fit/C_W. Reporting s as
        # well as C_fit is what makes his published values directly comparable.
        rows.append(dict(k=s.k, N=s.N, X=s.X, logX=s.logX, C_W=CW,
                         C_fit=r.C_fit, SE=r.SE, delta=r.C_fit - CW,
                         s_wolf=r.C_fit / CW, SE_s=r.SE / CW,
                         s_env_min=env.min() / CW, s_env_max=env.max() / CW,
                         g_max=int(r.g_max), n_points=int(r.n_points),
                         env_min=env.min(), env_max=env.max(),
                         env_spread=env.max() - env.min()))
    return pd.DataFrame(rows)


def grid_summary(grid: pd.DataFrame) -> pd.DataFrame:
    """Per specification: C_fit at the first and last scale, net change, monotone?"""
    out = []
    for (gmin, c, est), grp in grid.groupby(["g_min", "c", "estimator"], sort=True):
        grp = grp.sort_values("k")
        v = grp.C_fit.to_numpy()
        out.append(dict(g_min=gmin, c=c, estimator=est,
                        C_first=v[0], C_last=v[-1], net=v[-1] - v[0],
                        monotone_dec=bool(np.all(np.diff(v) < 0)),
                        monotone_inc=bool(np.all(np.diff(v) > 0))))
    return pd.DataFrame(out)


# ---------------------------------------------------------------- literature
# Wolf (arXiv:1102.0481) fits a(x)exp(-s(x)u) to log T_d(x) against
# u = g*pi(x)/x over the "linear portions" of the plot, and reports s falling
# from 1.187 to 1.136 as x runs from 2^28 to 2^48. Since s = C_fit/C_W, those
# values are directly comparable with this census over the overlapping range.
WOLF = {"s_lo": 1.187, "x_lo": 2 ** 28, "s_hi": 1.136, "x_hi": 2 ** 48,
        "cite": "Wolf, arXiv:1102.0481"}


def wolf_comparison(grid: pd.DataFrame, prim: pd.DataFrame) -> dict:
    """Compare the published slope with this census at the overlapping scale.

    Returns the scale nearest the start of Wolf's range, our value there, the
    envelope the specification grid generates, and the specification whose
    value is closest to his.
    """
    i = int(np.argmin(np.abs(np.log(prim.X.to_numpy()) - math.log(WOLF["x_lo"]))))
    row = prim.iloc[i]
    sub = grid[grid.k == row.k].copy()
    sub["s"] = sub.C_fit / row.C_W
    j = int((sub.s - WOLF["s_lo"]).abs().idxmin())
    best = sub.loc[j]
    return dict(X=int(row.X), logX=float(row.logX),
                s_primary=float(row.s_wolf),
                s_min=float(sub.s.min()), s_max=float(sub.s.max()),
                inside=bool(sub.s.min() <= WOLF["s_lo"] <= sub.s.max()),
                best_c=int(best.c), best_gmin=int(best.g_min),
                best_est=str(best.estimator), best_s=float(best.s),
                published=WOLF["s_lo"])
