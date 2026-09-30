"""
fgp.local -- window-free local diagnostics on the differenced intervals.

Differencing consecutive cumulative snapshots gives the exact gap distribution
on I_k = (X_{k-1}, X_k].  On each I_k we compute
    m1     mean gap
    E2     m2 / (2 m1^2)      (= 1 for an exponential law)
    E3     m3 / (6 m1^3)      (= 1 for an exponential law)
    D_left sup_g |S(g) - exp(-g/m1)|,   S the empirical survival function
    D_mid  the same with the midpoint convention exp(-(g+1)/m1)
    W2     Cramer-von Mises distance to the fitted exponential (added: a
           second discrepancy that weights the whole distribution, not just
           the point of maximal deviation)
plus multinomial bootstrap standard errors, consecutive-step z-scores
(independent under the multinomial model because the intervals are disjoint),
and a log-log power-law fit of each deviation against log X_eff.
"""

import math

import numpy as np
import pandas as pd
from scipy.special import expi

from .core import Scale, log_x_eff


def li(x: float) -> float:
    return expi(math.log(x))


def _stats(gs: np.ndarray, p: np.ndarray):
    """Moments and discrepancies. gs must be the FULL even grid up to the largest
    occupied gap, with p = 0 at unoccupied gaps: the survival function is a step
    function, so the supremum defining D can be attained at a gap that no prime
    carries, and restricting the sup to occupied gaps understates it."""
    m1 = float((gs * p).sum())
    m2 = float((gs ** 2 * p).sum())
    m3 = float((gs ** 3 * p).sum())
    S = 1.0 - np.cumsum(p)
    D0 = float(np.max(np.abs(S - np.exp(-gs / m1))))
    D1 = float(np.max(np.abs(S - np.exp(-(gs + 1.0) / m1))))
    # Cramer-von Mises: sum over cells of (F_emp - F_exp)^2 * p (discrete form)
    F = np.cumsum(p)
    Fe = 1.0 - np.exp(-(gs + 1.0) / m1)
    W2 = float(((F - Fe) ** 2 * p).sum())
    return m1, m2, m3, m2 / (2 * m1 ** 2), m3 / (6 * m1 ** 3), D0, D1, W2


def local_table(scales: list[Scale], B: int = 2000, seed: int = 20260803) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for prev, cur in zip(scales, scales[1:]):
        occ = np.union1d(prev.g, cur.g)
        gmax = int(occ.max())
        allg = np.arange(2, gmax + 2, 2)          # full even grid, holes included
        loc = np.array([cur.count_of(int(g)) - prev.count_of(int(g)) for g in allg],
                       dtype=np.int64)
        if (loc < 0).any():
            raise ValueError("differenced count is negative: censuses are not nested")
        gs, cs = allg.astype(float), loc.astype(float)
        n = int(cs.sum())
        p = cs / n
        m1, m2, m3, E2, E3, D0, D1, W2 = _stats(gs, p)

        bs = np.empty((B, 8))
        for b in range(B):
            bs[b] = _stats(gs, rng.multinomial(n, p) / n)
        se = bs.std(axis=0, ddof=1)

        lxe = log_x_eff(prev.X, cur.X)
        lxe_density = (cur.X - prev.X) / (li(cur.X) - li(prev.X))
        rows.append(dict(k=cur.k, X_lo=prev.X, X_hi=cur.X, n_gaps=n,
                         m1=m1, logX_eff=lxe, logX_eff_density=lxe_density,
                         ratio=lxe / m1,
                         E2=E2, SE_E2=se[3], E3=E3, SE_E3=se[4],
                         D_left=D0, SE_D_left=se[5], D_mid=D1, SE_D_mid=se[6],
                         W2=W2, SE_W2=se[7]))
    return pd.DataFrame(rows)


def step_tests(loc: pd.DataFrame) -> pd.DataFrame:
    """z-scores of consecutive steps, positive = toward the exponential law."""
    out = []
    specs = [("E2", "SE_E2", +1), ("E3", "SE_E3", +1),
             ("D_left", "SE_D_left", -1), ("D_mid", "SE_D_mid", -1), ("W2", "SE_W2", -1)]
    for name, sname, sign in specs:
        v, s = loc[name].to_numpy(), loc[sname].to_numpy()
        d = sign * np.diff(v)
        z = d / np.sqrt(s[1:] ** 2 + s[:-1] ** 2)
        end = sign * (v[-1] - v[0])
        end_se = math.sqrt(s[-1] ** 2 + s[0] ** 2)
        out.append(dict(statistic=name, n_steps=len(z), n_toward=int((z > 0).sum()),
                        min_z=float(z.min()), median_z=float(np.median(z)),
                        end_change=float(sign * end), end_SE=end_se,
                        end_z=float(end / end_se), z_steps=z))
    return pd.DataFrame(out)


def power_law_fits(loc: pd.DataFrame) -> pd.DataFrame:
    """Fit |stat - target| = A (log X_eff)^(-beta); report beta and extrapolations."""
    L = loc["logX_eff"].to_numpy()
    out = []
    for name, target in (("E2", 1.0), ("E3", 1.0), ("D_left", 0.0), ("D_mid", 0.0)):
        dev = np.abs(loc[name].to_numpy() - target)
        A = np.vstack([np.ones_like(L), np.log(L)]).T
        coef, *_ = np.linalg.lstsq(A, np.log(dev), rcond=None)
        resid = np.log(dev) - A @ coef
        s2 = float(resid @ resid) / (len(L) - 2)
        se = math.sqrt(s2 * np.linalg.inv(A.T @ A)[1, 1])
        pred = lambda LL: math.exp(coef[0] + coef[1] * math.log(LL))
        out.append(dict(statistic=name, dev_last=float(dev[-1]), beta=-coef[1],
                        SE_beta=se, logA=coef[0], at_50=pred(50), at_200=pred(200)))
    return pd.DataFrame(out)
