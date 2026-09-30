"""
fgp.spatial -- the relative spatial mean R(g;N), the crossover g*(X), and the
parameter-free density model.

R(g;N) = mu(g;N) / mu(N), mu(g;N) the mean of the primes in layer g and mu(N)
the mean of all classified primes p_2..p_N.

Crossover.  g*(X) is where R crosses 1.  Because R(g) is nearly flat near the
crossing and single layers carry arithmetic fluctuations, g* is located by a
linear regression of R on g over the band |R - 1| < b and solving for R = 1.
The band b is an analyst choice; we report b in {0.02, 0.03, 0.05} together
with plain linear interpolation between the two even gaps bracketing R = 1.

Density model.  Primes have density 1/log t.  The conditional probability that
a prime near t carries outgoing gap g is taken as
    'log'  :  (C2 H(g) / log t) exp(-g / log t)          (paper's form)
    'pi'   :  (C2 H(g) pi(t)/t) exp(-g pi(t) / t)          (Wolf's form)
In both cases C2 H(g) cancels in a ratio of means, so
    R_model(g;X) = [ int t w dt / int w dt ] / [ int t/log t dt / int 1/log t dt ],
w_g(t) = rate(t)^2 exp(-g rate(t)) / ... (see below).  Substituting t = e^{L-v}
keeps the quadrature stable for large L.
"""

import math

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expi

from .core import Scale

MIN_COUNT = 30       # layers with fewer primes are not used for R
BANDS = (0.02, 0.03, 0.05)          # bands in R (sensitivity)
XWINDOWS = ((0.7, 1.3), (0.6, 1.4), (0.5, 1.5))   # windows in x = g/log X
PRIMARY_XWINDOW = (0.6, 1.4)


def spatial_table(scales: list[Scale]) -> pd.DataFrame:
    """Layer means and the relative spatial mean R, with standard errors.

    R = mu_g / mu_all and mu_all contains mu_g, so the two are correlated.
    Writing w = n_g/N and splitting the primes into layer g and the rest,
    mu_all = w mu_g + (1-w) mu_r, and the delta method gives

        Var(R) = (1-w)^2 [ mu_r^2 sigma_g^2/n_g + mu_g^2 sigma_r^2/n_r ] / mu_all^4,

    which is computed exactly from the stored first and second moments. The
    standard errors are those of the sampling model in which layer membership is
    random; the primes themselves are deterministic, so they quantify the
    arithmetic fluctuation of a layer, not a probability about the primes.
    """
    rows = []
    for s in scales:
        tot_c = int(s.count.sum())
        tot_s = int(sum(int(v) for v in s.psum))       # exact, in Python integers
        tot_q = int(sum(int(v) for v in s.psqs)) if s.psqs is not None else None
        mu_all = tot_s / tot_c
        for i, (g, c, ps) in enumerate(zip(s.g, s.count, s.psum)):
            c = int(c)
            if c < MIN_COUNT:
                continue
            mu_g = int(ps) / c
            se = float("nan")
            if tot_q is not None and c >= 2:
                n_r = tot_c - c
                s1_r = tot_s - int(ps)
                s2_r = tot_q - int(s.psqs[i])
                if n_r >= 2:
                    mu_r = s1_r / n_r
                    var_r = (n_r * s2_r - s1_r * s1_r) / (n_r * (n_r - 1))
                    var_g = s.sd_of_index(i) ** 2
                    w = c / tot_c
                    v = (1 - w) ** 2 * (mu_r ** 2 * var_g / c
                                        + mu_g ** 2 * var_r / n_r) / mu_all ** 4
                    se = math.sqrt(v) if v > 0 else 0.0
            rows.append(dict(k=s.k, N=s.N, X=s.X, logX=s.logX, g=int(g),
                             x=g / s.logX, count=c, mu_g=mu_g, mu_all=mu_all,
                             R=mu_g / mu_all, SE_R=se))
    return pd.DataFrame(rows)


def _gstar_band(g, R, b):
    m = np.abs(R - 1) < b
    if m.sum() < 4:
        return float("nan")
    b1, b0 = np.polyfit(g[m], R[m], 1)
    return (1 - b0) / b1


def _gstar_interp(g, R):
    below = np.where(R < 1)[0]
    above = np.where(R >= 1)[0]
    if len(below) == 0 or len(above) == 0:
        return float("nan")
    i = below[-1]
    j = i + 1
    if j >= len(g):
        return float("nan")
    return g[i] + (1 - R[i]) * (g[j] - g[i]) / (R[j] - R[i])


def _gstar_xwindow(g, R, L, lo, hi):
    m = (g / L >= lo) & (g / L <= hi)
    if m.sum() < 4:
        return float("nan")
    b1, b0 = np.polyfit(g[m], R[m], 1)
    return (1 - b0) / b1


def _gstar_wls(g, R, se, L, lo, hi):
    """Crossover by weighted least squares, with a delta-method standard error.

    Fits R = b0 + b1 g with weights 1/SE^2 over the window, then solves R = 1.
    Var(g*) follows from the covariance of (b0, b1) by the delta method with
    dg*/db0 = -1/b1 and dg*/db1 = -(1-b0)/b1^2.
    """
    m = (g / L >= lo) & (g / L <= hi) & np.isfinite(se) & (se > 0)
    if m.sum() < 4:
        return float("nan"), float("nan")
    x, y, w = g[m], R[m], 1.0 / se[m] ** 2
    A = np.column_stack([np.ones_like(x), x])
    ATA = A.T @ (w[:, None] * A)
    b = np.linalg.solve(ATA, A.T @ (w * y))
    cov = np.linalg.inv(ATA)
    resid = y - A @ b
    chi2_red = float((w * resid ** 2).sum()) / (m.sum() - 2)
    cov = cov * max(chi2_red, 1.0)      # inflate if the linear form is imperfect
    b0, b1 = b
    gstar = (1 - b0) / b1
    d0, d1 = -1.0 / b1, -(1 - b0) / b1 ** 2
    var = d0 ** 2 * cov[0, 0] + d1 ** 2 * cov[1, 1] + 2 * d0 * d1 * cov[0, 1]
    return gstar, (math.sqrt(var) if var > 0 else float("nan"))


def crossover_table(spat: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for k, grp in spat.groupby("k", sort=True):
        grp = grp.sort_values("g")
        g, R = grp.g.to_numpy(float), grp.R.to_numpy()
        L = float(grp.logX.iloc[0])
        row = dict(k=k, X=int(grp.X.iloc[0]), logX=L)
        for lo, hi in XWINDOWS:
            row[f"gstar_x{int(lo*10)}{int(hi*10)}"] = _gstar_xwindow(g, R, L, lo, hi)
        for b in BANDS:
            row[f"gstar_b{int(b*100):02d}"] = _gstar_band(g, R, b)
        row["gstar_interp"] = _gstar_interp(g, R)
        lo, hi = PRIMARY_XWINDOW
        se = grp.SE_R.to_numpy()
        gw, sw = _gstar_wls(g, R, se, L, lo, hi)
        row["gstar_wls"] = gw
        row["SE_gstar"] = sw
        row["gstar"] = gw if np.isfinite(gw) else row[f"gstar_x{int(lo*10)}{int(hi*10)}"]
        row["ratio"] = row["gstar"] / row["logX"]
        row["deficit"] = row["logX"] - row["gstar"]
        rows.append(row)
    return pd.DataFrame(rows)


# ----------------------------------------------------------------- density model

def _li(x):
    return expi(math.log(x))


def _rate_u(u, variant):
    """Gap rate as a function of u = log t, without ever forming t."""
    if variant == "log":
        return 1.0 / u
    # pi(t)/t with pi(t) ~ li(t).  For u < 700, exp/expi are safe; beyond that use the
    # standard asymptotic expansion li(t)/t ~ (1/u)(1 + 1!/u + 2!/u^2 + ...).
    if u < 700.0:
        t = math.exp(u)
        return _li(t) / t
    return (1.0 + 1.0 / u + 2.0 / u ** 2 + 6.0 / u ** 3 + 24.0 / u ** 4) / u


def R_model(g: float, L: float, variant: str = "log") -> float:
    """Model relative spatial mean at log X = L for gap g."""
    # All four integrands carry a weight exp(-s v) with s >= 1, so the interval
    # beyond v = 50 contributes less than e^{-50} relative and truncating it keeps
    # the adaptive quadrature from under-sampling the bulk near v = 0 when L is large.
    hi = min(L - math.log(2.0), 50.0)

    def w(v, s):
        u = L - v                      # u = log t
        r = _rate_u(u, variant)
        # layer-g density is (1/log t) * r * exp(-g r); the factor e^{L-v} from
        # dt and the e^{L} common scale are handled by the e^{-s v} weights.
        return math.exp(-s * v) * r * math.exp(-g * r) / u

    A = quad(w, 0, hi, args=(2.0,), limit=300)[0]
    B = quad(w, 0, hi, args=(1.0,), limit=300)[0]
    C = quad(lambda v: math.exp(-2 * v) / (L - v), 0, hi, limit=300)[0]
    D = quad(lambda v: math.exp(-1 * v) / (L - v), 0, hi, limit=300)[0]
    return (A / B) / (C / D)


def g_star_model(L: float, variant: str = "log") -> float:
    return brentq(lambda g: R_model(g, L, variant) - 1.0, 1e-6, 5 * L, xtol=1e-10)


def model_asymptotics(variant: str = "log") -> pd.DataFrame:
    rows = []
    for e in (4, 6, 8, 10, 12, 16, 20, 30, 50, 100):
        L = e * math.log(10)
        g = g_star_model(L, variant)
        rows.append(dict(log10X=e, logX=L, gstar=g, ratio=g / L, deficit=L - g))
    df = pd.DataFrame(rows)
    # Richardson on the last two points: deficit = c0 + c1/L
    (L1, d1), (L2, d2) = df[["logX", "deficit"]].to_numpy()[-2:]
    c1 = (d1 - d2) / (1 / L1 - 1 / L2)
    df.attrs["richardson_c0"] = d2 - c1 / L2
    return df


def model_vs_data(spat: pd.DataFrame, ks: list[int], gmax_factor: float = 3.0,
                  min_count: int = 2000) -> pd.DataFrame:
    rows = []
    for k in ks:
        grp = spat[(spat.k == k) & (spat["count"] >= min_count)]
        L = float(grp.logX.iloc[0])
        grp = grp[grp.g <= gmax_factor * L]
        for _, r in grp.iterrows():
            rows.append(dict(k=k, X=int(r.X), logX=L, g=int(r.g), R=r.R,
                             R_log=R_model(r.g, L, "log"), R_pi=R_model(r.g, L, "pi")))
    return pd.DataFrame(rows)


def amplitude_ratio(spat: pd.DataFrame, lo: float = 0.3, hi: float = 1.5,
                    min_count: int = 2000) -> pd.DataFrame:
    """Slope of the collapsed profile log X (R-1) against x, census vs model.

    The census and the model are both fitted over the same window in the scaling
    variable x = g/log X, so the ratio of slopes measures how far the actual gap
    law is from the exponential envelope the model assumes. It is reported per
    scale, because quoting a residual in R makes the agreement look better than
    it is: R-1 is itself O(1/log X).
    """
    rows = []
    for k, g in spat.groupby("k", sort=True):
        L = float(g.logX.iloc[0])
        m = (g.x >= lo) & (g.x <= hi) & (g["count"] >= min_count)
        if int(m.sum()) < 5:
            continue
        x = g.x[m].to_numpy()
        sd = np.polyfit(x, L * (g.R[m].to_numpy() - 1), 1)[0]
        sm = np.polyfit(x, [L * (R_model(xx * L, L, "log") - 1) for xx in x], 1)[0]
        sp = np.polyfit(x, [L * (R_model(xx * L, L, "pi") - 1) for xx in x], 1)[0]
        rows.append(dict(k=int(k), X=int(g.X.iloc[0]), logX=L, slope_data=sd,
                         slope_log=sm, slope_pi=sp,
                         ratio_log=sd / sm, ratio_pi=sd / sp))
    return pd.DataFrame(rows)
