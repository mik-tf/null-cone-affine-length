import NullCone.Core
import NullCone.Existence

/-!
# Theorem 1 in the paper's own variables `(alpha, beta)`

`P = r_+ = lam0 / alpha`, `M = r_- = lam0 / beta`.  The two equations of Theorem 1 are exactly the
hypotheses `G*(lam0) = 1` and `int G* = Lambda` of the core lemma.
-/

open MeasureTheory intervalIntegral Set

namespace NullCone

noncomputable section

/-- `int_0^L G* = (P ln (P/(P-L)) - M ln ((M+L)/M)) / (P+M)` (partial fractions). -/
lemma integral_Gstar {L P M : ℝ} (hL : 0 < L) (hM : 0 < M) (hP : L < P) :
    ∫ l in (0:ℝ)..L, l / ((P - l) * (M + l))
      = (P * Real.log (P / (P - L)) - M * Real.log ((M + L) / M)) / (P + M) := by
  have hP0 : 0 < P := by linarith
  let H : ℝ → ℝ := fun x => (-(P * Real.log (P - x)) - M * Real.log (M + x)) / (P + M)
  have hH : ∀ x ∈ Icc (0:ℝ) L, HasDerivAt H (x / ((P - x) * (M + x))) x := by
    intro x hx
    obtain ⟨hx0, hx1⟩ := hx
    have h1 : HasDerivAt (fun y => M + y) 1 x := by simpa using (hasDerivAt_id x).const_add M
    have h2 : HasDerivAt (fun y => P - y) (-1) x := by simpa using (hasDerivAt_id x).const_sub P
    have l1 := (h1.log (by linarith)).const_mul M
    have l2 := (h2.log (by linarith)).const_mul P
    have := ((l2.neg).sub l1).div_const (P + M)
    convert this using 1
    have : P - x ≠ 0 := by linarith
    have : M + x ≠ 0 := by linarith
    field_simp
    ring
  have hcont : ContinuousOn (fun l => l / ((P - l) * (M + l))) (Icc 0 L) := by
    apply ContinuousOn.div continuousOn_id (by fun_prop)
    intro x hx
    have := Qf_pos (L := L) (P := P) hM hP hx
    exact this.ne'
  rw [integral_eq_sub_of_hasDerivAt (f := H) (a := 0) (b := L)
    (fun x hx => by rw [uIcc_of_le hL.le] at hx; exact hH x hx)
    (hcont.intervalIntegrable_of_Icc hL.le)]
  simp only [H]
  have hl1 : Real.log (P / (P - L)) = Real.log P - Real.log (P - L) :=
    Real.log_div hP0.ne' (by linarith)
  have hl2 : Real.log ((M + L) / M) = Real.log (M + L) - Real.log M :=
    Real.log_div (by linarith) hM.ne'
  rw [hl1, hl2]
  simp
  ring

/-- **Theorem 1, per-generator bound, in the paper's variables.**  Given `(alpha, beta)` solving
the two equations of Theorem 1 (with `lam0 = tau/eps`, `Lambda = ln (1/eps)`), every admissible `G`
(`C^1` on `[0, lam0]`, `G(lam0) = 1`, `int G = Lambda`, `G' >= G^2`, `G(0) >= 0`) has
`ell(G) <= D_max` and in fact `ell(G) <= D_max - (lam0 / sqrt (alpha beta)) G(0)`. -/
theorem theorem1_bound
    (L Lam α β : ℝ) (hL : 0 < L) (hα0 : 0 < α) (hα1 : α < 1) (hβ : 0 < β)
    (heq1 : L = α * β / ((1 - α) * (1 + β))) (heq2 : Lam = Lamf α β)
    (G G' : ℝ → ℝ)
    (hGc : ContinuousOn G (Icc 0 L))
    (hGd : ∀ x ∈ Ioo 0 L, HasDerivAt G (G' x) x)
    (hG'c : ContinuousOn G' (Icc 0 L))
    (h1 : G L = 1) (h2 : ∫ l in (0:ℝ)..L, G l = Lam)
    (h3 : ∀ l ∈ Icc 0 L, G l ^ 2 ≤ G' l) :
    ∫ l in (0:ℝ)..L, 2 * Real.sqrt (G' l - G l ^ 2)
      ≤ (2 * Real.sqrt (α * β) / (α + β)) * Real.log ((1 + β) / (1 - α))
        - (L / Real.sqrt (α * β)) * G 0 := by
  have h1a : 0 < 1 - α := by linarith
  have hαβ : 0 < α * β := mul_pos hα0 hβ
  set P := L / α with hPdef
  set M := L / β with hMdef
  have hM : 0 < M := by positivity
  have hLP : L < P := by
    rw [hPdef, lt_div_iff₀ hα0]; nlinarith
  have hPL : P - L = L * (1 - α) / α := by rw [hPdef]; field_simp
  have hML : M + L = L * (1 + β) / β := by rw [hMdef]; field_simp
  have hLM : L / M = β := by rw [hMdef]; field_simp
  have hLP' : L / P = α := by rw [hPdef]; field_simp
  have hs1 : L / ((P - L) * (M + L)) = 1 := by
    rw [hPL, hML]
    have : (1 - α) * (1 + β) ≠ 0 := by positivity
    field_simp
    rw [heq1] at *
    field_simp
    ring_nf
    sorry
  sorry

end

end NullCone
