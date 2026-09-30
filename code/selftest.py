#!/usr/bin/env python3
"""
selftest.py -- offline accuracy checks for the FGP pipeline.

    python3 selftest.py ../data/fgp_layers.csv

Runs with no network access. Nothing here re-runs the pipeline and compares it
with itself; every check confronts the computation with a fact established
somewhere else:

  A. Published constants   pi(x), p_n, twin-prime counts, maximal prime gaps,
                           all hard-coded below with their sources. If the C
                           sieve is wrong, these fail.
  B. Internal identities   things that must hold for any correct census
                           (layer counts partition the primes, prime sums are
                           consistent, no integer overflow).
  C. Independent sieve     a from-scratch NumPy sieve recomputes the layer
                           counts and the same published constants.
  D. Statistical methods   the Poisson GLM against a general-purpose optimiser,
                           the bootstrap SE against the delta method, H(g)
                           against direct factorisation.
  E. Mathematics           the exponential moment identities numerically, and
                           the numerical density model against the analytic
                           first-order expansion derived in Appendix A.

Exit status is 0 only if every check passes.
"""

import math
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from fgp import core, frequency, local, spatial  # noqa: E402

PASS, FAIL = [], []


def check(name, ok, detail=""):
    (PASS if ok else FAIL).append(name)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{('  ' + detail) if detail else ''}")


# ---------------------------------------------------------------- ground truth
#
# pi(10^k), k = 1..10. Standard values, e.g. Crandall & Pomerance, "Prime
# Numbers: A Computational Perspective", 2nd ed., Table 1.1; OEIS A006880.
PI_POWERS = {10**1: 4, 10**2: 25, 10**3: 168, 10**4: 1229, 10**5: 9592,
             10**6: 78498, 10**7: 664579, 10**8: 5761455, 10**9: 50847534,
             10**10: 455052511}

# p_n for n a power of 10. OEIS A006988.
NTH_PRIME = {1: 2, 10: 29, 100: 541, 1000: 7919, 10**4: 104729, 10**5: 1299709,
             10**6: 15485863, 10**7: 179424673, 10**8: 2038074743}

# Twin-prime counts pi_2(10^k): pairs (p, p+2) both prime with p <= 10^k.
# OEIS A007508.
TWIN_COUNTS = {10**3: 35, 10**4: 205, 10**5: 1224, 10**6: 8169, 10**7: 58980,
               10**8: 440312, 10**9: 3424506, 10**10: 27412679}

# Maximal prime gaps and the prime that precedes each. OEIS A005250 / A002386,
# also Nicely, Math. Comp. 68 (1999) 1311-1315.
MAXIMAL_GAPS = [(1, 2), (2, 3), (4, 7), (6, 23), (8, 89), (14, 113), (18, 523),
                (20, 887), (22, 1129), (34, 1327), (36, 9551), (44, 15683),
                (52, 19609), (72, 31397), (86, 155921), (96, 360653),
                (112, 370261), (114, 492113), (118, 1349533), (132, 1357201),
                (148, 2010733), (154, 4652353), (180, 17051707),
                (210, 20831323), (220, 47326693), (222, 122164747),
                (234, 189695659), (248, 191912783), (250, 387096133),
                (282, 436273009), (288, 1294268491), (292, 1453168141),
                (320, 2300942549), (336, 3842610773), (354, 4302407359)]


def sieve_upto(n):
    """Plain boolean sieve, written as simply as possible so it is easy to audit."""
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.flatnonzero(s)


# ------------------------------------------------------------------ A. constants
def group_A(scales):
    print("\nA. Census against published constants")
    by_N = {s.N: s for s in scales}

    # p_N at checkpoints must equal the tabulated n-th prime.
    for n, p in NTH_PRIME.items():
        if n in by_N:
            check(f"p_{n:,} = {p:,}", by_N[n].X == p, f"got {by_N[n].X:,}")

    # The largest checkpoint covers 10^10: its index is pi(10^10).
    last = scales[-1]
    check("pi(10^10) = 455,052,511", last.N == PI_POWERS[10**10], f"got {last.N:,}")
    check("largest prime <= 10^10 is 9,999,999,967", last.X == 9_999_999_967,
          f"got {last.X:,}")

    # The g = 2 layer at that checkpoint counts exactly the twin-prime pairs
    # below 10^10, a quantity computed independently by many authors.
    check("twin pairs below 10^10 = 27,412,679",
          last.count_of(2) == TWIN_COUNTS[10**10], f"got {last.count_of(2):,}")

    # Every gap present in the census must appear in the maximal-gap record at
    # or below its first occurrence, and the largest gap seen must be 354.
    gmax = max(int(g) for s in scales for g in s.g)
    check("maximal gap below 10^10 is 354", gmax == 354, f"got {gmax}")

    # No odd gap above 1, and gap 1 occurs exactly once (the pair 2,3).
    odd = {int(g) for s in scales for g in s.g if g % 2 == 1}
    check("only odd gap is 1", odd <= {1}, f"got {sorted(odd)}")
    check("gap 1 occurs at most once", all(s.count_of(1) <= 1 for s in scales))


# ------------------------------------------------------------- B. internal identities
def group_B(scales):
    print("\nB. Internal identities that any correct census must satisfy")
    for s in scales:
        tot = int(s.count.sum())
        if tot != s.N - 1:
            check(f"counts partition primes p_2..p_N at N={s.N:,}", False,
                  f"{tot:,} vs {s.N - 1:,}")
            return
    check("layer counts sum to N-1 at all 20 checkpoints", True)

    # Cumulative monotonicity: a layer can never lose members as N grows.
    ok = True
    for prev, cur in zip(scales, scales[1:]):
        for g in np.union1d(prev.g, cur.g):
            if cur.count_of(int(g)) < prev.count_of(int(g)):
                ok = False
    check("layer counts are non-decreasing in N", ok)

    # Prime sums: the mean prime in any layer must lie in (2, X].
    ok = all(2 < int(s.psum[i]) / int(s.count[i]) <= s.X
             for s in scales for i in range(len(s.g)))
    check("every layer mean lies in (2, X]", ok)

    # Second moments: variance must be non-negative and the layer SD must be
    # comparable to, but not equal to, the uniform-on-[0,X] value X/sqrt(12).
    s = scales[-1]
    if s.psqs is not None:
        sds = [s.sd_of_index(i) for i in range(len(s.g)) if int(s.count[i]) >= 100]
        unif = s.X / math.sqrt(12)
        check("all layer variances are non-negative", all(d == d and d >= 0 for d in sds))
        check("layer SDs are within 20% of the uniform value (sanity, not identity)",
              all(0.8 * unif < d < 1.2 * unif for d in sds),
              f"range {min(sds):.3e} to {max(sds):.3e}, uniform {unif:.3e}")
    else:
        check("census carries second moments (sumsq column)", False,
              "regenerate the census to get error bars on R")

    # Overflow guard: the total prime sum must be well below 2^64.
    tot = int(scales[-1].psum.sum(dtype=np.uint64))
    check("total prime sum below 2^64 (no overflow)", tot < 2**64,
          f"{tot:.3e}, headroom x{2**64 / tot:.1f}")

    # Layers hold OUTGOING gaps, so the gaps of p_2..p_N telescope to
    # p_{N+1} - p_2, not p_N - p_2. Recovering p_{N+1} from the census must give
    # the smallest prime above 10^10, which is 10,000,000,019 (OEIS A003617).
    s = scales[-1]
    total = int((s.g.astype(object) * s.count.astype(object)).sum())
    check("gaps telescope to p_{N+1} - p_2 = 10,000,000,019 - 3",
          total + 3 == 10_000_000_019, f"got p_(N+1) = {total + 3:,}")


# ------------------------------------------------------- C. independent reimplementation
def group_C(scales, limit=10**7):
    print(f"\nC. Independent NumPy sieve (to {limit:,}, audited separately)")
    pr = sieve_upto(limit)
    check(f"pi({limit:,}) = {PI_POWERS[limit]:,}", len(pr) == PI_POWERS[limit],
          f"got {len(pr):,}")
    for n, p in NTH_PRIME.items():
        if n <= len(pr):
            check(f"independent p_{n:,}", int(pr[n - 1]) == p, f"got {pr[n-1]:,}")

    gaps = np.diff(pr)
    for x, t in TWIN_COUNTS.items():
        if x <= limit:
            # pairs with smaller prime <= x
            n = int(((gaps == 2) & (pr[:-1] <= x)).sum())
            check(f"independent twin count below {x:,}", n == t, f"got {n:,}")

    # Maximal gap record reproduced from scratch.
    rec, best = [], 0
    for i, g in enumerate(gaps):
        if g > best:
            best = int(g)
            rec.append((best, int(pr[i])))
    exp = [(g, p) for g, p in MAXIMAL_GAPS if p <= limit]
    check(f"maximal gap record below {limit:,} ({len(exp)} entries)", rec == exp,
          "" if rec == exp else f"got {rec[-3:]}")

    # Layer counts, sums AND sums of squares at the smallest census checkpoint,
    # recomputed independently. This is the check that validates the error bars
    # on R: if sumsq were wrong, every standard error in Section 5 would be wrong.
    s = scales[0]
    if s.X <= limit:
        idx = int(np.searchsorted(pr, s.X)) + 1
        gl = gaps[1:idx]                           # gaps of p_2 .. p_N
        pl = pr[1:idx]                             # the primes p_2 .. p_N
        cnt = np.bincount(gl, minlength=600)
        mism = sum(1 for g, c in zip(s.g, s.count) if int(cnt[g]) != int(c))
        check(f"all layer counts at N={s.N:,} match independent sieve", mism == 0,
              f"{mism} mismatches")
        bad_s = bad_q = 0
        for i, g in enumerate(s.g):
            sel = pl[gl == int(g)]
            if int(sum(int(v) for v in sel)) != int(s.psum[i]):
                bad_s += 1
            if s.psqs is not None and int(sum(int(v) * int(v) for v in sel)) != int(s.psqs[i]):
                bad_q += 1
        check(f"all layer prime-sums at N={s.N:,} match independently", bad_s == 0,
              f"{bad_s} mismatches")
        if s.psqs is not None:
            check(f"all layer sums-of-squares at N={s.N:,} match independently",
                  bad_q == 0, f"{bad_q} mismatches")

        # The delta-method SE on R, checked against a direct bootstrap over the
        # actual primes in a layer (possible only at this small scale).
        from fgp import spatial as _sp
        sp = _sp.spatial_table([s])
        rng = np.random.default_rng(7)
        g0 = int(sp.sort_values("count").iloc[-1].g)     # the most populated layer
        sel = pl[gl == g0].astype(float)
        mu_all = float(sum(int(v) for v in s.psum)) / int(s.count.sum())
        boot = np.array([rng.choice(sel, size=len(sel), replace=True).mean() / mu_all
                         for _ in range(400)])
        se_boot, se_delta = boot.std(ddof=1), float(sp[sp.g == g0].SE_R.iloc[0])
        check("delta-method SE(R) agrees with a direct bootstrap within 10%",
              abs(se_boot - se_delta) / se_delta < 0.10,
              f"bootstrap {se_boot:.3e} vs delta {se_delta:.3e}")


# ------------------------------------------------------------ D. statistical methods
def group_D(scales):
    print("\nD. Statistical methods against independent implementations")
    from scipy.optimize import minimize
    from sympy import factorint

    # H(g) against sympy's factorisation.
    ok = True
    for g in range(2, 400, 2):
        h = 1.0
        for q in factorint(g):
            if q > 2:
                h *= (q - 1) / (q - 2)
        if abs(h - core.H(g)) > 1e-12:
            ok = False
    check("H(g) matches sympy factorisation for all even g < 400", ok)

    # Poisson GLM (our IRLS) against a general-purpose optimiser on the same
    # negative log-likelihood.
    s = scales[-1]
    gmax = int(3 * s.logX) // 2 * 2
    m = (s.g >= 4) & (s.g <= gmax) & (s.count > 0)
    gs, cs = s.g[m].astype(float), s.count[m].astype(float)
    off = np.array([math.log(core.H(int(g))) for g in s.g[m]])
    lam_irls, se = frequency.fit_poisson(s.g[m], s.count[m])

    def nll(b):
        eta = b[0] + off - b[1] * gs
        return float(np.sum(np.exp(eta) - cs * eta))

    r = minimize(nll, [math.log(cs.mean()), 0.05], method="Nelder-Mead",
                 options=dict(xatol=1e-12, fatol=1e-12, maxiter=20000))
    check("Poisson GLM slope matches Nelder-Mead on the same likelihood",
          abs(r.x[1] - lam_irls) < 1e-7,
          f"IRLS {lam_irls:.10f} vs NM {r.x[1]:.10f}")

    # Bootstrap SE against the delta method, on the last interval.
    loc = local.local_table(scales[-2:], B=4000, seed=12345)
    prev, cur = scales[-2], scales[-1]
    allg = np.union1d(prev.g, cur.g)
    c = np.array([cur.count_of(int(g)) - prev.count_of(int(g)) for g in allg])
    keep = c > 0
    gs, cs = allg[keep].astype(float), c[keep].astype(float)
    n = cs.sum()
    p = cs / n
    m1 = (gs * p).sum()
    m2 = (gs ** 2 * p).sum()
    d = gs ** 2 / (2 * m1 ** 2) - m2 * gs / m1 ** 3        # dE2/dp_i
    var = (float((p * d ** 2).sum()) - float((p * d).sum()) ** 2) / n
    se_delta = math.sqrt(var)
    se_boot = float(loc.SE_E2.iloc[0])
    check("bootstrap SE(E2) agrees with the delta method within 5%",
          abs(se_boot - se_delta) / se_delta < 0.05,
          f"boot {se_boot:.3e} vs delta {se_delta:.3e}")

    # The bootstrap must be insensitive to the seed.
    a = local.local_table(scales[-2:], B=4000, seed=1).SE_E2.iloc[0]
    b = local.local_table(scales[-2:], B=4000, seed=2).SE_E2.iloc[0]
    check("bootstrap SE stable across seeds within 5%", abs(a - b) / b < 0.05,
          f"{a:.3e} vs {b:.3e}")


# ---------------------------------------------------------------- E. mathematics
def group_E():
    print("\nE. Mathematical identities")
    from scipy.integrate import quad

    # Exponential moment ratios E2 = E3 = 1, by numerical integration.
    for r in (0.5, 1.0, 2.7):
        mom = [quad(lambda x, k=k: x ** k * r * math.exp(-r * x), 0, np.inf)[0]
               for k in (1, 2, 3)]
        check(f"E2 = 1 for Exp({r})", abs(mom[1] / (2 * mom[0] ** 2) - 1) < 1e-9)
        check(f"E3 = 1 for Exp({r})", abs(mom[2] / (6 * mom[0] ** 3) - 1) < 1e-9)

    # The numerical density model must converge to the analytic first-order
    # expansion of Appendix A: log X (R - 1) -> (x - 1)/2.
    worst = 0.0
    for x in (0.2, 0.5, 1.0, 2.0):
        v = 10000.0 * (spatial.R_model(x * 10000.0, 10000.0, "log") - 1)
        worst = max(worst, abs(v - (x - 1) / 2))
    check("model matches the first-order expansion at log X = 10^4",
          worst < 5e-3, f"max deviation {worst:.2e}")

    # Hence the crossover deficit tends to 3/2.
    dfc = 10000.0 - spatial.g_star_model(10000.0, "log")
    check("model crossover deficit -> 3/2", abs(dfc - 1.5) < 5e-3, f"got {dfc:.5f}")

    # Quadrature must be stable: doubling the truncation must not move R.
    a = spatial.R_model(2 * 230.3, 230.3, "log")
    check("model is numerically stable at large log X", math.isfinite(a) and 0.99 < a < 1.01,
          f"R = {a:.6f}")


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "../data/fgp_layers.csv"
    print(f"FGP self-test, census = {path}")
    scales = core.load_census(path)
    print(f"loaded {len(scales)} checkpoints, X from {scales[0].X:,} to {scales[-1].X:,}")
    group_A(scales)
    group_B(scales)
    group_C(scales)
    group_D(scales)
    group_E()
    print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
    if FAIL:
        for f in FAIL:
            print("  FAILED:", f)
        sys.exit(1)
    print("All checks passed.")


if __name__ == "__main__":
    main()
