/-
FGP.lean -- elementary statements from "Forward-Gap Layers of Consecutive Primes"
formalised in Lean 4 with Mathlib.

STATUS: every Mathlib name used below was checked to exist, by name and
signature, against Mathlib at revision 18dc857edd2da9c0d63672b94fbaf48a9f4810b5
(Lean v4.35.0-rc3).  One name was wrong and has been corrected: the constant
pulled out of an integral is `integral_const_mul`, not `integral_mul_left`.
The proof TERMS have NOT been machine-checked end to end: the prebuilt Mathlib
cache was unreachable from the machine where this was prepared, and compiling
Mathlib from source there was not practical.  Run `lake exe cache get && lake
build` (see the project README) before citing Appendix B as verified.  The four
results are elementary; the remaining risk is in the proof scripts, not the
statements.

Contents
  1. Forward-gap layers are pairwise disjoint and every index 2 <= n <= N lies
     in exactly one layer.
  2. Boundary lemma: for every prime q <= X the next prime is <= 2X
     (Bertrand's postulate), so sieving [2, 2X] classifies every prime <= X.
  3. Telescoping: the forward gaps of p_2, ..., p_N sum to p_{N+1} - p_2.
  4. For the exponential law, E[Z^2] = 2 (E Z)^2 and E[Z^3] = 6 (E Z)^3,
     i.e. the moment ratios E_2 and E_3 of the paper equal 1.
-/

import Mathlib

open Finset

namespace FGP

/-! ### 1. Layers -/

/-- The `n`-th prime, 1-indexed: `p 1 = 2`, `p 2 = 3`, ...  Mathlib's `Nat.nth` is 0-indexed. -/
noncomputable def p (n : ℕ) : ℕ := Nat.nth Nat.Prime (n - 1)

/-- The forward gap `d n = p (n+1) - p n`. -/
noncomputable def d (n : ℕ) : ℕ := p (n + 1) - p n

/-- The forward-gap layer, recorded by prime index: indices `n` with `2 <= n <= N` and `d n = g`. -/
noncomputable def layer (g N : ℕ) : Finset ℕ :=
  (Finset.Icc 2 N).filter (fun n => d n = g)

theorem mem_layer {g N n : ℕ} : n ∈ layer g N ↔ (2 ≤ n ∧ n ≤ N) ∧ d n = g := by
  simp only [layer, Finset.mem_filter, Finset.mem_Icc]

/-- Layers with different gaps are disjoint. -/
theorem layer_disjoint {g g' N : ℕ} (h : g ≠ g') : Disjoint (layer g N) (layer g' N) := by
  unfold layer
  rw [Finset.disjoint_filter]
  intro n _ h1 h2
  exact h (h1.symm.trans h2)

/-- Every index `2 <= n <= N` lies in exactly one layer, namely the one indexed by its own gap. -/
theorem exists_unique_layer {N n : ℕ} (hn : 2 ≤ n ∧ n ≤ N) : ∃! g, n ∈ layer g N := by
  refine ⟨d n, ?_, ?_⟩
  · exact mem_layer.mpr ⟨hn, rfl⟩
  · intro g hg
    exact (mem_layer.mp hg).2.symm

/-! ### 2. Boundary lemma -/

/-- For every prime `q <= X` there is a prime `r` with `q < r <= 2X`, and the least such `r`
(the next prime after `q`) is therefore at most `2X`.  Consequently a sieve of `[2, 2X]`
records the outgoing gap of every prime up to `X`.  Uses Bertrand's postulate
`Nat.exists_prime_lt_and_le_two_mul`, verified present with signature
`(n : ℕ) (hn0 : n ≠ 0) : ∃ p, Nat.Prime p ∧ n < p ∧ p ≤ 2 * n`. -/
theorem next_prime_le_two_mul {X q : ℕ} (hq : q.Prime) (hqX : q ≤ X) :
    ∃ r, r.Prime ∧ q < r ∧ r ≤ 2 * X ∧ ∀ s, s.Prime → q < s → r ≤ s := by
  classical
  obtain ⟨r₀, hr₀, hqr₀, hr₀2⟩ := Nat.exists_prime_lt_and_le_two_mul q hq.ne_zero
  have hex : ∃ r, r.Prime ∧ q < r := ⟨r₀, hr₀, hqr₀⟩
  refine ⟨Nat.find hex, (Nat.find_spec hex).1, (Nat.find_spec hex).2, ?_, ?_⟩
  · exact (Nat.find_min' hex ⟨hr₀, hqr₀⟩).trans (hr₀2.trans (Nat.mul_le_mul_left 2 hqX))
  · intro s hs hqs
    exact Nat.find_min' hex ⟨hs, hqs⟩

/-! ### 3. Telescoping -/

/-- The forward gaps of `p 2, ..., p N` sum to `p (N+1) - p 2` (stated over `ℤ` to avoid
truncated subtraction).  With `N - 1` gaps this is the identity behind reading the Wolf
benchmark `π(X) log X / X` as a method-of-moments estimator. -/
theorem gap_telescope {N : ℕ} (hN : 1 ≤ N) :
    ∑ i in Finset.range (N - 1), ((p (i + 3) : ℤ) - p (i + 2)) = (p (N + 1) : ℤ) - p 2 := by
  -- Finset.sum_range_sub is the to_additive image of Finset.prod_range_div.
  have h := Finset.sum_range_sub (fun i => (p (i + 2) : ℤ)) (N - 1)
  have h2 : N - 1 + 2 = N + 1 := by omega
  simpa [h2, add_assoc] using h

/-! ### 4. Exponential moments -/

/-- The `k`-th raw moment of the exponential law with rate `r`, written as an integral. -/
noncomputable def expMoment (r : ℝ) (k : ℕ) : ℝ :=
  ∫ x in Set.Ioi (0 : ℝ), x ^ (k : ℝ) * (r * Real.exp (-r * x))

/-- `E[Z^k] = k! / r^k` for `Z ~ Exp(r)`.  Uses `integral_rpow_mul_exp_neg_mul_Ioi`,
`Real.Gamma_nat_eq_factorial` and `MeasureTheory.integral_const_mul`. -/
theorem expMoment_eq {r : ℝ} (hr : 0 < r) (k : ℕ) :
    expMoment r k = (k.factorial : ℝ) / r ^ k := by
  unfold expMoment
  have hk : (-1 : ℝ) < k := by
    have : (0 : ℝ) ≤ k := Nat.cast_nonneg k
    linarith
  have hI := integral_rpow_mul_exp_neg_mul_Ioi hk hr
  -- pull the constant r out of the integrand
  have hint :
      (∫ x in Set.Ioi (0 : ℝ), x ^ (k : ℝ) * (r * Real.exp (-r * x)))
        = r * ∫ x in Set.Ioi (0 : ℝ), x ^ (k : ℝ) * Real.exp (-r * x) := by
    rw [← MeasureTheory.integral_const_mul]
    congr 1
    funext x
    ring
  rw [hint, hI]
  have hG : Real.Gamma ((k : ℝ) + 1) = (k.factorial : ℝ) := by
    exact_mod_cast Real.Gamma_nat_eq_factorial k
  rw [hG]
  have hpow : (1 / r) ^ ((k : ℝ) + 1) = 1 / r ^ (k + 1) := by
    rw [show ((k : ℝ) + 1) = ((k + 1 : ℕ) : ℝ) by push_cast; ring, Real.rpow_natCast]
    rw [one_div, inv_pow, one_div]
  rw [hpow]
  field_simp
  ring

/-- `E_2 = m_2 / (2 m_1^2) = 1` for the exponential law. -/
theorem exp_E2 {r : ℝ} (hr : 0 < r) :
    expMoment r 2 / (2 * (expMoment r 1) ^ 2) = 1 := by
  rw [expMoment_eq hr 2, expMoment_eq hr 1]
  simp only [Nat.factorial, Nat.cast_ofNat, Nat.cast_one]
  field_simp
  ring

/-- `E_3 = m_3 / (6 m_1^3) = 1` for the exponential law. -/
theorem exp_E3 {r : ℝ} (hr : 0 < r) :
    expMoment r 3 / (6 * (expMoment r 1) ^ 3) = 1 := by
  rw [expMoment_eq hr 3, expMoment_eq hr 1]
  simp only [Nat.factorial, Nat.cast_ofNat, Nat.cast_one]
  field_simp
  ring

end FGP
