"""
fgp.verify -- an independent reimplementation of the layer census.

This is deliberately different from fgp_sieve.c: NumPy, one byte per integer
(not one bit per odd integer), segments of 10^8, gaps taken with np.diff rather
than with a running 'prev'.  It recomputes |L(g;N)| at every checkpoint with
X <= limit and compares with the C census exactly.  A mismatch anywhere is a
bug in one of the two programs.
"""

import math
import time

import numpy as np
import pandas as pd

from .core import Scale


def _base_primes(n: int) -> np.ndarray:
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.flatnonzero(s)


def independent_census(scales: list[Scale], limit: int, seg: int = 100_000_000,
                       gmax: int = 512, margin: int = 1000) -> pd.DataFrame:
    """Return a table of (k, g, count_C, count_py, match) for checkpoints with X <= limit.

    The sieve runs to ``limit + margin`` rather than to ``limit``. A prime's
    layer is fixed by the NEXT prime, so a run stopping exactly at ``limit``
    cannot classify its own last prime and therefore cannot report the
    checkpoint at N = pi(limit) at all: this is the boundary convention of
    Section 2.2, and the verifier is subject to it exactly as the census is.
    With limit = 10^10 the next prime is 10,000,000,019, so any margin above 19
    suffices; 1000 covers the maximal gap comfortably at this size.
    """
    checks = {s.N: s for s in scales if s.X <= limit}
    top = limit + margin
    base = _base_primes(int(math.isqrt(top)) + 1)
    cnt = np.zeros(gmax, dtype=np.int64)
    carry = None            # last prime of the previous segment
    n_prev = 0              # index of `carry` in the prime sequence (1-based)
    results = {}
    t0 = time.time()
    lo = 2
    while lo <= top:
        hi = min(lo + seg, top + 1)
        s = np.ones(hi - lo, dtype=bool)
        for p in base:
            if p * p >= hi:
                break
            start = max(p * p, ((lo + p - 1) // p) * p)
            s[start - lo::p] = False
        if lo == 2:
            pass
        else:
            s[:max(0, 2 - lo)] = False
        primes = np.flatnonzero(s) + lo
        if carry is not None:
            primes = np.concatenate(([carry], primes))
        gaps = np.diff(primes)
        # prime primes[i] has index n_prev + i + (1 if carry is None else 0)
        # we need layer contributions for indices n >= 2 (skip p_1 = 2)
        idx0 = n_prev if carry is not None else 1   # index of primes[0]
        idx = idx0 + np.arange(len(gaps))           # index of each gap's prime
        ok = idx >= 2
        if ok.any() and int(gaps[ok].max()) >= gmax:
            j = int(np.flatnonzero(ok & (gaps >= gmax))[0])
            raise ValueError(
                f"gap {int(gaps[j])} at prime {int(primes[j]):,} exceeds gmax="
                f"{gmax}; raise gmax rather than truncating the histogram")
        cnt += np.bincount(gaps[ok], minlength=gmax)[:gmax]
        # checkpoints inside this segment: N with idx0 <= N <= idx0+len(gaps)-1
        for N, sc in checks.items():
            if N in results:
                continue
            if idx0 <= N <= idx0 + len(gaps) - 1:
                # counts up to and including gap of prime index N: subtract later gaps
                later = idx > N
                snap = cnt - np.bincount(gaps[ok & later], minlength=gmax)[:gmax]
                results[N] = snap.copy()
        carry = primes[-1]
        n_prev = idx0 + len(primes) - 1
        lo = hi
    missing = sorted(set(checks) - set(results))
    if missing:
        raise RuntimeError(
            f"checkpoints {missing} were never reached: the sieve stopped at "
            f"{top:,}, which is not far enough past X to supply the outgoing "
            f"gap of the last classified prime. Increase `margin`.")
    rows = []
    for N, sc in checks.items():
        snap = results[N]
        for g, c in zip(sc.g, sc.count):
            rows.append(dict(k=sc.k, N=N, X=sc.X, g=int(g), count_C=int(c),
                             count_py=int(snap[g]) if g < gmax else -1))
        # also any g present in python but absent in C
        for g in np.flatnonzero(snap):
            if sc.count_of(int(g)) == 0:
                rows.append(dict(k=sc.k, N=N, X=sc.X, g=int(g), count_C=0,
                                 count_py=int(snap[g])))
    df = pd.DataFrame(rows).sort_values(["k", "g"]).reset_index(drop=True)
    df["match"] = df.count_C == df.count_py
    df.attrs["seconds"] = time.time() - t0
    return df
