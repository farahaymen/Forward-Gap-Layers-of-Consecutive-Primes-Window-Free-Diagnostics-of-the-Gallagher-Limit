"""
fgp.core -- load the layer census and provide the arithmetic factor H(g).

The census file (produced by fgp_sieve.c) has columns k,N,X,g,count,sum:
    k      checkpoint index (0..19)
    N      prime index at the checkpoint, X = p_N
    X      the N-th prime
    g      even gap
    count  |L(g;N)|, the number of primes p_n (2 <= n <= N) with d_n = g
    sum    sum of those primes (exact, 64-bit)

Released under CC0 / public domain.
"""

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

C2 = 1.3203236316  # twin-prime constant, 2 * prod_{q>2} (1 - 1/(q-1)^2)


@dataclass
class Scale:
    k: int
    N: int
    X: int
    g: np.ndarray       # even gaps with nonzero count, ascending
    count: np.ndarray   # |L(g;N)|
    psum: np.ndarray    # sum of primes in L(g;N)
    psqs: np.ndarray | None = None   # sum of squares (None for legacy censuses)

    @property
    def logX(self) -> float:
        return math.log(self.X)

    @property
    def ngaps(self) -> int:
        return int(self.count.sum())

    def count_of(self, g: int) -> int:
        i = np.searchsorted(self.g, g)
        return int(self.count[i]) if i < len(self.g) and self.g[i] == g else 0

    def sum_of(self, g: int) -> int:
        i = np.searchsorted(self.g, g)
        return int(self.psum[i]) if i < len(self.g) and self.g[i] == g else 0

    def sd_of_index(self, i: int) -> float:
        """Sample SD of the primes in layer i, from the stored moments.

        Needs the sumsq column; returns nan for a legacy census. The moments are
        exact integers, so the subtraction is done in Python integers to avoid
        the catastrophic cancellation a float64 evaluation would suffer
        (sum^2/n and sumsq agree to within ~1e-2 relative at these scales).
        """
        if self.psqs is None:
            return float("nan")
        n = int(self.count[i])
        if n < 2:
            return float("nan")
        s1, s2 = int(self.psum[i]), int(self.psqs[i])
        var = (n * s2 - s1 * s1) / (n * (n - 1))   # exact numerator, then divide
        return math.sqrt(var) if var > 0 else 0.0


def load_census(path: str) -> list[Scale]:
    """Load a census. The sumsq column exceeds 64 bits, so it is read as exact
    Python integers via dtype=object rather than as a numpy integer type."""
    dt = {"k": int, "N": int, "X": int, "g": int, "count": int, "sum": "uint64"}
    head = pd.read_csv(path, nrows=0)
    has_sq = "sumsq" in head.columns
    if has_sq:
        dt["sumsq"] = object
    df = pd.read_csv(path, dtype=dt)
    if has_sq:
        df["sumsq"] = df["sumsq"].map(int)
    scales = []
    for k, grp in df.groupby("k", sort=True):
        grp = grp.sort_values("g")
        scales.append(Scale(
            k=int(k), N=int(grp["N"].iloc[0]), X=int(grp["X"].iloc[0]),
            g=grp["g"].to_numpy(dtype=np.int64),
            count=grp["count"].to_numpy(dtype=np.int64),
            psum=grp["sum"].to_numpy(dtype=np.uint64),
            psqs=grp["sumsq"].to_numpy(dtype=object) if has_sq else None,
        ))
    return scales


def H(g: int) -> float:
    """Hardy--Littlewood arithmetic factor prod_{q | g, q > 2 prime} (q-1)/(q-2).

    The factor 2 must be divided out first. Leaving it in makes the trial-division
    remainder even, so the closing branch applies (m-1)/(m-2) to a composite even
    number as though it were an odd prime; it also lets the loop bound q*q <= m
    terminate before small odd factors are found (for g = 6 the loop never runs).
    """
    if g <= 0:
        return 1.0
    h, m = 1.0, g
    while m % 2 == 0:
        m //= 2
    q = 3
    while q * q <= m:
        if m % q == 0:
            h *= (q - 1) / (q - 2)
            while m % q == 0:
                m //= q
        q += 2
    if m > 1:
        h *= (m - 1) / (m - 2)
    return h


def log_x_eff(Xlo: int, Xhi: int) -> float:
    """Mean of log t over (Xlo, Xhi] weighted by the prime density 1/log t.

    int_{Xlo}^{Xhi} (log t)(1/log t) dt / int_{Xlo}^{Xhi} dt/log t
    = (Xhi - Xlo) / li-difference.  The closed form used in the paper,
    (Xhi(log Xhi - 1) - Xlo(log Xlo - 1)) / (Xhi - Xlo), is the mean of log t
    with respect to Lebesgue measure and agrees with the density-weighted mean
    to O(1/log X); both are tabulated by run_all so the difference is visible.
    """
    return (Xhi * (math.log(Xhi) - 1) - Xlo * (math.log(Xlo) - 1)) / (Xhi - Xlo)
