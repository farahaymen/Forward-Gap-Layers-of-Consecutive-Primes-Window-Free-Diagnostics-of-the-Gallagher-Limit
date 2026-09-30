#!/usr/bin/env python3
"""
make_census.py -- produce the layer census WITHOUT a C compiler.

This is a drop-in replacement for fgp_sieve.c: same checkpoint schedule, same
CSV columns, byte-identical output.  It is slower than the C sieve (a few
minutes to 10^10 instead of ~15 s) but needs nothing but Python and NumPy, so
the whole pipeline runs on Windows or anywhere else without a build toolchain.

    python make_census.py --limit 10000001000 --pi 455052511 --out ../data/fgp_layers.csv

Smaller runs for a quick check (the analysis works on any of these, the numbers
just correspond to a shorter range):

    python make_census.py --limit 100000000  --out ../data/fgp_layers_1e8.csv
    python make_census.py --limit 1000000000 --out ../data/fgp_layers_1e9.csv

Note on --limit: pass a limit slightly ABOVE the largest X you want classified,
because a prime's layer is determined by the NEXT prime.  The production run
uses 10^10 + 10^3.  If --pi is omitted it is computed on the fly, which costs
nothing extra but means the checkpoint schedule is built after the first pass;
passing it explicitly (as the C sieve requires) reproduces the shipped file.

Released under CC0 / public domain.
"""

import argparse
import math
import sys
import time

import numpy as np

GMAX = 512
SEG = 1 << 26          # 64 Mi integers per segment; ~64 MB as a bool array


def build_schedule(pi_max: int) -> list[int]:
    """The same 20-point geometric schedule in prime index N as fgp_sieve.c."""
    K, lo, hi = 20, 2.0e5, float(pi_max)
    r = (hi / lo) ** (1.0 / (K - 1))
    ck, v = [], lo
    for _ in range(K):
        ck.append(int(v + 0.5))
        v *= r
    ck[K - 1] = pi_max
    for forced in (1_000_000, 9_000_000, 100_000_000):
        best, bd = 0, 1e300
        for k in range(K - 1):
            d = abs(math.log(ck[k]) - math.log(forced))
            if d < bd:
                bd, best = d, k
        ck[best] = forced
    return sorted(ck)


def base_primes(n: int) -> np.ndarray:
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.flatnonzero(s).astype(np.int64)


def segments(limit: int, seg: int):
    """Yield arrays of primes in increasing order, segment by segment."""
    base = base_primes(int(math.isqrt(limit)) + 1)
    lo = 2
    while lo <= limit:
        hi = min(lo + seg, limit + 1)
        s = np.ones(hi - lo, dtype=bool)
        for p in base:
            if p * p >= hi:
                break
            start = max(p * p, ((lo + p - 1) // p) * p)
            s[start - lo::p] = False
        if lo <= 1:
            s[:2 - lo] = False
        yield np.flatnonzero(s).astype(np.int64) + lo
        lo = hi


def count_primes(limit: int, seg: int, verbose: bool) -> int:
    n, t0 = 0, time.time()
    for primes in segments(limit, seg):
        n += len(primes)
    if verbose:
        print(f"  pi({limit}) = {n:,}  ({time.time() - t0:.1f}s)", file=sys.stderr)
    return n


def census_exact(limit: int, pi_max: int, seg: int, out, verbose: bool):
    """Stream the primes, accumulating layer counts and EXACT layer prime-sums.

    The sums are accumulated in Python integers rather than floats: they reach
    about 2e18 at X = 10^10, past the 2^53 where float64 stops being exact, and
    the spatial statistic R(g;N) needs them exact.
    """
    ck = build_schedule(pi_max)
    cnt = np.zeros(GMAX, dtype=object)
    ssum = np.zeros(GMAX, dtype=object)
    ssqs = np.zeros(GMAX, dtype=object)
    cnt[:] = 0
    ssum[:] = 0
    ssqs[:] = 0
    snap = {}
    todo = list(ck)
    carry, idx_carry = None, 0
    t0 = time.time()

    for primes in segments(limit, seg):
        if carry is not None:
            primes = np.concatenate(([carry], primes))
            idx0 = idx_carry
        else:
            idx0 = 1
        gaps = np.diff(primes)
        if len(gaps) == 0:
            carry, idx_carry = primes[-1], idx0 + len(primes) - 1
            continue
        idx = idx0 + np.arange(len(gaps), dtype=np.int64)
        too_big = (idx >= 2) & (gaps >= GMAX)
        if too_big.any():
            j = int(np.flatnonzero(too_big)[0])
            raise SystemExit(f"FATAL: gap {int(gaps[j])} at prime {int(primes[j]):,} "
                             f"exceeds GMAX={GMAX}; raise GMAX and rerun")
        ok = idx >= 2
        g_ok, p_ok, i_ok = gaps[ok], primes[:-1][ok], idx[ok]

        def accumulate(gs, ps, c, s, q):
            """Add the (gap, prime) pairs into count/sum/sumsq arrays, exactly."""
            order = np.argsort(gs, kind="stable")
            gs, ps = gs[order], ps[order]
            uniq, starts = np.unique(gs, return_index=True)
            ends = np.append(starts[1:], len(gs))
            for g, a, b in zip(uniq, starts, ends):
                c[g] += int(b - a)
                # int(ps[a:b].sum()) would overflow int64 only above 9.2e18;
                # segment sums stay well below that, so this is exact.
                s[g] += int(ps[a:b].sum(dtype=np.uint64))
                # p^2 reaches 1e20, past uint64, so square in Python integers.
                q[g] += sum(int(v) * int(v) for v in ps[a:b])
            return c, s, q

        while todo and idx0 <= todo[0] <= idx0 + len(gaps) - 1:
            N = todo.pop(0)
            upto = i_ok <= N
            c2, s2, q2 = cnt.copy(), ssum.copy(), ssqs.copy()
            accumulate(g_ok[upto], p_ok[upto], c2, s2, q2)
            snap[N] = (c2, s2, q2, int(primes[N - idx0]))
            if verbose:
                print(f"  checkpoint N={N:,} X={snap[N][3]:,} ({time.time() - t0:.1f}s)",
                      file=sys.stderr)

        accumulate(g_ok, p_ok, cnt, ssum, ssqs)
        carry, idx_carry = primes[-1], idx0 + len(primes) - 1

    out.write("k,N,X,g,count,sum,sumsq\n")
    rows = 0
    for k, N in enumerate(ck):
        if N not in snap:
            print(f"warning: checkpoint N={N:,} never reached (limit too small)", file=sys.stderr)
            continue
        c, s, q, X = snap[N]
        for g in range(GMAX):
            if c[g]:
                out.write(f"{k},{N},{X},{g},{c[g]},{s[g]},{q[g]}\n")
                rows += 1
    if verbose:
        print(f"  wrote {rows} rows in {time.time() - t0:.1f}s", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description="Layer census without a C compiler.")
    ap.add_argument("--limit", type=int, default=10_000_001_000,
                    help="sieve limit; use slightly above the largest X you want (default 10^10+10^3)")
    ap.add_argument("--pi", type=int, default=0,
                    help="pi(limit), if known; omitted means one extra counting pass")
    ap.add_argument("--out", default="../data/fgp_layers.csv")
    ap.add_argument("--segment", type=int, default=SEG, help="segment size in integers")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    verbose = not a.quiet

    pi_max = a.pi
    if not pi_max:
        if verbose:
            print("counting primes to build the checkpoint schedule...", file=sys.stderr)
        pi_max = count_primes(a.limit, a.segment, verbose)
        # the schedule tops out at pi(X) for the largest X classified, which is
        # pi(limit) minus the primes in the overshoot; the C sieve is passed
        # pi(10^10) = 455052511 for limit = 10^10 + 10^3.

    if verbose:
        print(f"census to {a.limit:,} with pi_max = {pi_max:,}", file=sys.stderr)
    with open(a.out, "w", newline="") as fh:
        census_exact(a.limit, pi_max, a.segment, fh, verbose)
    if verbose:
        print(f"wrote {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
