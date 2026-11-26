import Mathlib

/-!
# Theorem 1: existence of `(alpha, beta)` (intermediate value theorem)

With `L = lam0 = tau/eps` and `Lam = ln (1/eps)`: the paper's condition `tau > 1 - eps` is
`Lam < ln (1 + L)`, and `eps < 1` is `0 < Lam`.
The curve is parametrised by `r = L (1 - alpha) / alpha` in `(0, 1)`, that is
`alpha = L / (L + r)`, `beta = r / (1 - r)`, so the first equation of the theorem holds identically.
-/

open Set Filter Topology

namespace NullCone

noncomputable section

/-- The right side of the second equation of Theorem 1, as a function of `(alpha, beta)`. -/
def Lamf (α β : ℝ) : ℝ := (β * Real.log (1 / (1 - α)) - α * Real.log (1 + β)) / (α + β)

/-- The same function along the curve, divided through by `beta`, which is continuous up to `r = 1`. -/
def fr (L r : ℝ) : ℝ :=
  (Real.log ((L + r) / r) + L / ((L + r) * r) * ((1 - r) * Real.log (1 - r))) /
    (L / ((L + r) * r) * (1 - r) + 1)

/-- The same function along the curve, written to be continuous up to `r = 0`. -/
def f0 (L r : ℝ) : ℝ :=
  ((r * Real.log (L + r) - r * Real.log r) / (1 - r) + L / (L + r) * Real.log (1 - r)) /
    (L / (L + r) + r / (1 - r))

lemma fr_eq_Lamf {L r : ℝ} (hL : 0 < L) (hr0 : 0 < r) (hr1 : r < 1) :
    fr L r = Lamf (L / (L + r)) (r / (1 - r)) := by
  have h1 : 0 < 1 - r := by linarith
  have h2 : 0 < L + r := by linarith
  have e1 : 1 - L / (L + r) = r / (L + r) := by field_simp; ring
  have e2 : 1 / (r / (L + r)) = (L + r) / r := by field_simp
  have e3 : 1 + r / (1 - r) = 1 / (1 - r) := by field_simp; ring
  unfold fr Lamf
  rw [e1, e2, e3, one_div, Real.log_inv]
  have hp : 0 < L / ((L + r) * r) * (1 - r) + 1 := by positivity
  have hp2 : 0 < L / (L + r) + r / (1 - r) := by positivity
  rw [div_eq_div_iff hp.ne' hp2.ne']
  field_simp
  ring

lemma fr_eq_f0 {L r : ℝ} (hL : 0 < L) (hr0 : 0 < r) (hr1 : r < 1) : fr L r = f0 L r := by
  have h1 : 0 < 1 - r := by linarith
  have h2 : 0 < L + r := by linarith
  have e : Real.log ((L + r) / r) = Real.log (L + r) - Real.log r :=
    Real.log_div h2.ne' hr0.ne'
  unfold fr f0
  rw [e]
  have hp : 0 < L / ((L + r) * r) * (1 - r) + 1 := by positivity
  have hp2 : 0 < L / (L + r) + r / (1 - r) := by positivity
  rw [div_eq_div_iff hp.ne' hp2.ne']
  field_simp
  ring

lemma continuousOn_fr {L : ℝ} (hL : 0 < L) : ContinuousOn (fr L) (Ioo 0 1) := by
  intro r hr
  obtain ⟨hr0, hr1⟩ := hr
  have h1 : 0 < 1 - r := by linarith
  have h2 : 0 < L + r := by linarith
  have hcm : Continuous fun x : ℝ => x * Real.log x := Real.continuous_mul_log
  have hc2 : ContinuousAt (fun x : ℝ => (1 - x) * Real.log (1 - x)) r :=
    (hcm.comp (by fun_prop : Continuous fun x : ℝ => 1 - x)).continuousAt
  have hk : ContinuousAt (fun x : ℝ => L / ((L + x) * x)) r := by
    apply ContinuousAt.div continuousAt_const (by fun_prop) (by positivity)
  have hlog : ContinuousAt (fun x : ℝ => Real.log ((L + x) / x)) r := by
    apply ContinuousAt.log (ContinuousAt.div (by fun_prop) continuousAt_id hr0.ne')
    positivity
  have hnum : ContinuousAt (fun x : ℝ => Real.log ((L + x) / x)
      + L / ((L + x) * x) * ((1 - x) * Real.log (1 - x))) r := hlog.add (hk.mul hc2)
  have hden : ContinuousAt (fun x : ℝ => L / ((L + x) * x) * (1 - x) + 1) r :=
    (hk.mul (by fun_prop)).add continuousAt_const
  have hp : L / ((L + r) * r) * (1 - r) + 1 ≠ 0 := by positivity
  exact (hnum.div hden hp).continuousWithinAt

lemma tendsto_fr_one {L : ℝ} (hL : 0 < L) :
    Tendsto (fr L) (𝓝[<] 1) (𝓝 (Real.log (1 + L))) := by
  have hcm : Continuous fun x : ℝ => x * Real.log x := Real.continuous_mul_log
  have hc2 : ContinuousAt (fun x : ℝ => (1 - x) * Real.log (1 - x)) 1 :=
    (hcm.comp (by fun_prop : Continuous fun x : ℝ => 1 - x)).continuousAt
  have hk : ContinuousAt (fun x : ℝ => L / ((L + x) * x)) 1 := by
    apply ContinuousAt.div continuousAt_const (by fun_prop) (by positivity)
  have hlog : ContinuousAt (fun x : ℝ => Real.log ((L + x) / x)) 1 := by
    apply ContinuousAt.log (ContinuousAt.div (by fun_prop) continuousAt_id one_ne_zero)
    positivity
  have hnum : ContinuousAt (fun x : ℝ => Real.log ((L + x) / x)
      + L / ((L + x) * x) * ((1 - x) * Real.log (1 - x))) 1 := hlog.add (hk.mul hc2)
  have hden : ContinuousAt (fun x : ℝ => L / ((L + x) * x) * (1 - x) + 1) 1 :=
    (hk.mul (by fun_prop)).add continuousAt_const
  have hp : L / ((L + 1) * 1) * (1 - 1) + 1 ≠ 0 := by simp
  have hc : ContinuousAt (fr L) 1 := hnum.div hden hp
  have := hc.tendsto
  have e : fr L 1 = Real.log (1 + L) := by
    simp [fr, add_comm]
  rw [e] at this
  exact this.mono_left nhdsWithin_le_nhds

lemma tendsto_fr_zero {L : ℝ} (hL : 0 < L) : Tendsto (fr L) (𝓝[>] 0) (𝓝 0) := by
  have hcm : Continuous fun x : ℝ => x * Real.log x := Real.continuous_mul_log
  have hlog1 : ContinuousAt (fun x : ℝ => Real.log (L + x)) 0 :=
    ContinuousAt.log (by fun_prop) (by simpa using hL.ne')
  have hlog2 : ContinuousAt (fun x : ℝ => Real.log (1 - x)) 0 :=
    ContinuousAt.log (by fun_prop) (by simp)
  have ha : ContinuousAt (fun x : ℝ => L / (L + x)) 0 :=
    ContinuousAt.div continuousAt_const (by fun_prop) (by simpa using hL.ne')
  have hb : ContinuousAt (fun x : ℝ => (x * Real.log (L + x) - x * Real.log x) / (1 - x)) 0 :=
    ContinuousAt.div
      ((continuousAt_id.mul hlog1).sub hcm.continuousAt) (by fun_prop) (by simp)
  have hnum : ContinuousAt (fun x : ℝ => (x * Real.log (L + x) - x * Real.log x) / (1 - x)
      + L / (L + x) * Real.log (1 - x)) 0 := hb.add (ha.mul hlog2)
  have hden : ContinuousAt (fun x : ℝ => L / (L + x) + x / (1 - x)) 0 :=
    ha.add (ContinuousAt.div continuousAt_id (by fun_prop) (by simp))
  have hp : L / (L + 0) + 0 / (1 - 0) ≠ 0 := by simp [hL.ne']
  have hc : ContinuousAt (f0 L) 0 := hnum.div hden hp
  have e0 : f0 L 0 = 0 := by simp [f0, hL.ne']
  have ht := hc.tendsto
  rw [e0] at ht
  have ht' := ht.mono_left (nhdsWithin_le_nhds (s := Ioi (0:ℝ)))
  refine ht'.congr' ?_
  filter_upwards [Ioo_mem_nhdsGT (zero_lt_one : (0:ℝ) < 1)] with r hr
  exact (fr_eq_f0 hL hr.1 hr.2).symm

/-- **Existence of `(alpha, beta)`** (paper, Theorem 1, "Existence of `(alpha, beta)`"). -/
theorem exists_alpha_beta {L Lam : ℝ} (hL : 0 < L) (hLam0 : 0 < Lam)
    (hLam : Lam < Real.log (1 + L)) :
    ∃ α β : ℝ, 0 < α ∧ α < 1 ∧ 0 < β ∧
      L = α * β / ((1 - α) * (1 + β)) ∧ Lam = Lamf α β := by
  have hpre : IsPreconnected (Ioo (0:ℝ) 1) := isPreconnected_Ioo
  have hmem : Lam ∈ Ioo (0:ℝ) (Real.log (1 + L)) := ⟨hLam0, hLam⟩
  obtain ⟨r, ⟨hr0, hr1⟩, hr⟩ := hpre.intermediate_value_Ioo
    (l₁ := 𝓝[>] (0:ℝ)) (l₂ := 𝓝[<] (1:ℝ))
    (by rw [le_principal_iff]; exact Ioo_mem_nhdsGT zero_lt_one)
    (by rw [le_principal_iff]; exact Ioo_mem_nhdsLT zero_lt_one)
    (continuousOn_fr hL) (tendsto_fr_zero hL) (tendsto_fr_one hL) hmem
  have h1 : 0 < 1 - r := by linarith
  have h2 : 0 < L + r := by linarith
  refine ⟨L / (L + r), r / (1 - r), by positivity, ?_, by positivity, ?_, ?_⟩
  · rw [div_lt_one h2]; linarith
  · have : 1 - L / (L + r) = r / (L + r) := by field_simp; ring
    rw [this]
    have : 1 + r / (1 - r) = 1 / (1 - r) := by field_simp; ring
    rw [this]
    field_simp
  · rw [← fr_eq_Lamf hL hr0 hr1]; exact hr.symm

end

end NullCone
