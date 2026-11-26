import Mathlib

/-!
# Theorem 1: the core inequality

Notation (paper -> Lean): `lam0 = L`, `r_+ = P`, `r_- = M`, `Lambda = Lam`,
`Q(l) = (P - l) (M + l)`, `G*(l) = l / Q(l)`, `c = 2 sqrt (P M)`.
-/

open MeasureTheory intervalIntegral Set

namespace NullCone

noncomputable section

/-- `Q(l) = (r_+ - l)(r_- + l)`. -/
def Qf (P M l : ℝ) : ℝ := (P - l) * (M + l)

/-- The candidate `G*(l) = l / Q(l)`. -/
def Gstar (P M l : ℝ) : ℝ := l / Qf P M l

/-- The pointwise tangent bound, linear in `(g, q)`: AM-GM on `u s` and `(u g - l)^2 >= 0`.
Here `u = Q(l)`, `m = sqrt (P M)`, `g = G`, `q = G'`. -/
lemma ptwise (m u l g q : ℝ) (hg : g ^ 2 ≤ q) :
    2 * m * u * Real.sqrt (q - g ^ 2) ≤ u ^ 2 * q - 2 * l * u * g + l ^ 2 + m ^ 2 := by
  have hs := Real.sq_sqrt (sub_nonneg.2 hg)
  set s := Real.sqrt (q - g ^ 2) with hsdef
  have hq : u ^ 2 * q = u ^ 2 * g ^ 2 + u ^ 2 * s ^ 2 := by
    have : q = g ^ 2 + s ^ 2 := by linarith
    rw [this]; ring
  nlinarith [sq_nonneg (u * g - l), sq_nonneg (u * s - m)]

lemma hasDerivAt_Qf (P M x : ℝ) : HasDerivAt (Qf P M) (P - M - 2 * x) x := by
  have h1 : HasDerivAt (fun y => P - y) (-1) x := by simpa using (hasDerivAt_id x).const_sub P
  have h2 : HasDerivAt (fun y => M + y) 1 x := by simpa using (hasDerivAt_id x).const_add M
  have h := h1.mul h2
  have e : Qf P M = fun y => (P - y) * (M + y) := rfl
  rw [e]
  convert h using 1
  ring

lemma hasDerivAt_Gstar (P M x : ℝ) (hx : Qf P M x ≠ 0) :
    HasDerivAt (Gstar P M) ((x ^ 2 + P * M) / (Qf P M x) ^ 2) x := by
  have h := (hasDerivAt_id x).div (hasDerivAt_Qf P M x) hx
  have e : Gstar P M = fun y => id y / Qf P M y := rfl
  rw [e]
  convert h using 1
  simp only [id]
  unfold Qf
  field_simp
  ring

lemma Qf_pos {L P M x : ℝ} (hM : 0 < M) (hP : L < P) (hx : x ∈ Icc 0 L) :
    0 < Qf P M x := by
  obtain ⟨h0, h1⟩ := hx
  unfold Qf
  exact mul_pos (by linarith) (by linarith)

lemma continuous_Qf (P M : ℝ) : Continuous (Qf P M) := by
  unfold Qf; fun_prop

lemma continuousOn_Gstar {L P M : ℝ} (hM : 0 < M) (hP : L < P) :
    ContinuousOn (Gstar P M) (Icc 0 L) := by
  unfold Gstar
  exact continuousOn_id.div (continuous_Qf P M).continuousOn
    (fun x hx => (Qf_pos hM hP hx).ne')

/-- The closed form of `int 2 m / Q`. -/
lemma integral_inv_Q {L P M : ℝ} (hL : 0 < L) (hM : 0 < M) (hP : L < P) (m : ℝ) :
    ∫ l in (0:ℝ)..L, 2 * m / Qf P M l
      = (2 * m / (P + M)) * Real.log ((1 + L / M) / (1 - L / P)) := by
  have hP0 : 0 < P := by linarith
  let K : ℝ → ℝ := fun x => (2 * m / (P + M)) * (Real.log (M + x) - Real.log (P - x))
  have hK : ∀ x ∈ Icc (0:ℝ) L, HasDerivAt K (2 * m / Qf P M x) x := by
    intro x hx
    obtain ⟨hx0, hx1⟩ := hx
    have h1 : HasDerivAt (fun y => M + y) 1 x := by simpa using (hasDerivAt_id x).const_add M
    have h2 : HasDerivAt (fun y => P - y) (-1) x := by simpa using (hasDerivAt_id x).const_sub P
    have l1 := h1.log (by linarith)
    have l2 := h2.log (by linarith)
    have := (l1.sub l2).const_mul (2 * m / (P + M))
    convert this using 1
    have hPx : P - x ≠ 0 := by linarith
    have hMx : M + x ≠ 0 := by linarith
    have hPM : P + M ≠ 0 := by linarith
    unfold Qf
    field_simp
    ring
  have hcont : ContinuousOn (fun l => 2 * m / Qf P M l) (Icc 0 L) :=
    continuousOn_const.div (continuous_Qf P M).continuousOn
      (fun x hx => (Qf_pos hM hP hx).ne')
  have := integral_eq_sub_of_hasDerivAt (f := K) (f' := fun l => 2 * m / Qf P M l)
    (a := 0) (b := L)
    (fun x hx => by rw [uIcc_of_le hL.le] at hx; exact hK x hx)
    (hcont.intervalIntegrable_of_Icc hL.le)
  rw [this]
  simp only [K]
  have e : (1 + L / M) / (1 - L / P) = ((M + L) * P) / ((P - L) * M) := by
    have : P - L ≠ 0 := by linarith
    field_simp
  rw [e, Real.log_div (by positivity) (by nlinarith),
    Real.log_mul (by linarith) hP0.ne', Real.log_mul (by linarith) hM.ne']
  simp
  ring

/-- **Core lemma, strong form.**  In the paper's notation this is
`ell(G) <= ell(G*) - (c/2) G(0)` with `ell(G*) = D_max`. -/
theorem core_strong
    (L P M Lam : ℝ) (hL : 0 < L) (hM : 0 < M) (hP : L < P)
    (hstar1 : L / ((P - L) * (M + L)) = 1)
    (hstar2 : ∫ l in (0:ℝ)..L, l / ((P - l) * (M + l)) = Lam)
    (G G' : ℝ → ℝ)
    (hGc : ContinuousOn G (Icc 0 L))
    (hGd : ∀ x ∈ Ioo 0 L, HasDerivAt G (G' x) x)
    (hG'c : ContinuousOn G' (Icc 0 L))
    (h1 : G L = 1) (h2 : ∫ l in (0:ℝ)..L, G l = Lam)
    (h3 : ∀ l ∈ Icc 0 L, G l ^ 2 ≤ G' l) :
    ∫ l in (0:ℝ)..L, 2 * Real.sqrt (G' l - G l ^ 2)
      ≤ (2 * Real.sqrt (P * M) / (P + M)) * Real.log ((1 + L / M) / (1 - L / P))
        - Real.sqrt (P * M) * G 0 := by
  have hP0 : 0 < P := by linarith
  have hPM : 0 < P * M := mul_pos hP0 hM
  set m := Real.sqrt (P * M) with hmdef
  have hm : 0 < m := Real.sqrt_pos.2 hPM
  have hm2 : m ^ 2 = P * M := Real.sq_sqrt hPM.le
  have hLe : uIcc (0:ℝ) L = Icc 0 L := uIcc_of_le hL.le
  have hQpos : ∀ x ∈ Icc (0:ℝ) L, 0 < Qf P M x := fun x hx => Qf_pos hM hP hx
  have hQc : Continuous (Qf P M) := continuous_Qf P M
  have hGsc := continuousOn_Gstar (L := L) hM hP
  -- eta = G - G*
  let η : ℝ → ℝ := fun x => G x - Gstar P M x
  have hηc : ContinuousOn η (Icc 0 L) := hGc.sub hGsc
  -- derivative of Gstar written with the same shape
  let Gd : ℝ → ℝ := fun x => (x ^ 2 + P * M) / (Qf P M x) ^ 2
  have hGdc : ContinuousOn Gd (Icc 0 L) := by
    apply ContinuousOn.div (by fun_prop) (hQc.pow 2).continuousOn
    intro x hx; exact pow_ne_zero 2 (hQpos x hx).ne'
  -- Phi = Q * eta and its derivative
  let Φ : ℝ → ℝ := fun x => Qf P M x * η x
  let Φd : ℝ → ℝ := fun x => (P - M - 2 * x) * η x + Qf P M x * (G' x - Gd x)
  have hΦc : ContinuousOn Φ (Icc 0 L) := hQc.continuousOn.mul hηc
  have hΦd : ∀ x ∈ Ioo (0:ℝ) L, HasDerivAt Φ (Φd x) x := by
    intro x hx
    have hxI : x ∈ Icc (0:ℝ) L := Ioo_subset_Icc_self hx
    have := (hasDerivAt_Qf P M x).mul
      ((hGd x hx).sub (hasDerivAt_Gstar P M x (hQpos x hxI).ne'))
    exact this
  have hΦdc : ContinuousOn Φd (Icc 0 L) := by
    have a1 : ContinuousOn (fun x => (P - M - 2 * x) * η x) (Icc 0 L) :=
      (by fun_prop : Continuous fun x : ℝ => P - M - 2 * x).continuousOn.mul hηc
    exact a1.add (hQc.continuousOn.mul (hG'c.sub hGdc))
  have hFTC : ∫ x in (0:ℝ)..L, Φd x = Φ L - Φ 0 :=
    integral_eq_sub_of_hasDeriv_right_of_le hL.le hΦc
      (fun x hx => (hΦd x hx).hasDerivWithinAt)
      (hΦdc.intervalIntegrable_of_Icc hL.le)
  have hΦL : Φ L = 0 := by
    simp only [Φ, η, Gstar]
    have : Qf P M L = (P - L) * (M + L) := rfl
    rw [this, h1, hstar1]; ring
  have hΦ0 : Φ 0 = P * M * G 0 := by
    simp [Φ, η, Gstar, Qf]
  -- integral of eta is zero
  have hGint : IntervalIntegrable G volume 0 L := hGc.intervalIntegrable_of_Icc hL.le
  have hGsint : IntervalIntegrable (Gstar P M) volume 0 L := hGsc.intervalIntegrable_of_Icc hL.le
  have hη0 : ∫ x in (0:ℝ)..L, η x = 0 := by
    simp only [η]
    rw [integral_sub hGint hGsint, h2]
    have : ∫ x in (0:ℝ)..L, Gstar P M x = Lam := by
      rw [← hstar2]; rfl
    rw [this]; ring
  -- the pointwise bound
  let R : ℝ → ℝ := fun x => 2 * m / Qf P M x + (1 / m) * (Φd x - (P - M) * η x)
  have hRc : ContinuousOn R (Icc 0 L) := by
    apply ContinuousOn.add
    · exact continuousOn_const.div hQc.continuousOn (fun x hx => (hQpos x hx).ne')
    · exact continuousOn_const.mul (hΦdc.sub (continuousOn_const.mul hηc))
  have hbound : ∀ x ∈ Icc (0:ℝ) L, 2 * Real.sqrt (G' x - G x ^ 2) ≤ R x := by
    intro x hx
    have hu := hQpos x hx
    have hp := ptwise m (Qf P M x) x (G x) (G' x) (h3 x hx)
    have hmu : 0 < m * Qf P M x := mul_pos hm hu
    have key : 2 * Real.sqrt (G' x - G x ^ 2)
        ≤ (Qf P M x ^ 2 * G' x - 2 * x * Qf P M x * G x + x ^ 2 + m ^ 2) / (m * Qf P M x) := by
      rw [le_div_iff₀ hmu]; nlinarith
    refine key.trans (le_of_eq ?_)
    simp only [R, Φd, η, Gd, Gstar]
    have hu0 : Qf P M x ≠ 0 := hu.ne'
    field_simp
    nlinarith [hm2]
  have hsqc : ContinuousOn (fun x => 2 * Real.sqrt (G' x - G x ^ 2)) (Icc 0 L) := by
    apply ContinuousOn.mul continuousOn_const
    exact Real.continuous_sqrt.comp_continuousOn (hG'c.sub (hGc.pow 2))
  have hmono := intervalIntegral.integral_mono_on (μ := volume) hL.le
    (hsqc.intervalIntegrable_of_Icc hL.le) (hRc.intervalIntegrable_of_Icc hL.le) hbound
  refine hmono.trans (le_of_eq ?_)
  -- evaluate the integral of R
  have hint1 : IntervalIntegrable (fun x => 2 * m / Qf P M x) volume 0 L :=
    (continuousOn_const.div hQc.continuousOn (fun x hx => (hQpos x hx).ne')).intervalIntegrable_of_Icc hL.le
  have hint2 : IntervalIntegrable (fun x => Φd x - (P - M) * η x) volume 0 L :=
    (hΦdc.sub (continuousOn_const.mul hηc)).intervalIntegrable_of_Icc hL.le
  have hint3 : IntervalIntegrable Φd volume 0 L := hΦdc.intervalIntegrable_of_Icc hL.le
  have hint4 : IntervalIntegrable (fun x => (P - M) * η x) volume 0 L :=
    (continuousOn_const.mul hηc).intervalIntegrable_of_Icc hL.le
  show ∫ x in (0:ℝ)..L, R x = _
  simp only [R]
  rw [intervalIntegral.integral_add hint1 (hint2.const_mul _),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_sub hint3 hint4,
    intervalIntegral.integral_const_mul, hη0, hFTC, hΦL, hΦ0, integral_inv_Q hL hM hP]
  have : (1 / m) * (P * M * G 0) = m * G 0 := by
    have h5 : P * M * G 0 = m * (m * G 0) := by rw [← hm2]; ring
    rw [h5, ← mul_assoc, one_div_mul_cancel hm.ne', one_mul]
  have e3 : (1 / m) * (0 - P * M * G 0 - (P - M) * 0) = - (m * G 0) := by
    have : (1 / m) * (0 - P * M * G 0 - (P - M) * 0) = -((1 / m) * (P * M * G 0)) := by ring
    rw [this]; congr 1
  rw [e3]; ring

/-- **Core lemma, as stated in the feasibility note.**  `ell(G) <= D_max`. -/
theorem core_lemma
    (L P M Lam : ℝ) (hL : 0 < L) (hM : 0 < M) (hP : L < P)
    (hstar1 : L / ((P - L) * (M + L)) = 1)
    (hstar2 : ∫ l in (0:ℝ)..L, l / ((P - l) * (M + l)) = Lam)
    (G G' : ℝ → ℝ)
    (hGc : ContinuousOn G (Icc 0 L))
    (hGd : ∀ x ∈ Ioo 0 L, HasDerivAt G (G' x) x)
    (hG'c : ContinuousOn G' (Icc 0 L))
    (h1 : G L = 1) (h2 : ∫ l in (0:ℝ)..L, G l = Lam)
    (h3 : ∀ l ∈ Icc 0 L, G l ^ 2 ≤ G' l) (h4 : 0 ≤ G 0) :
    ∫ l in (0:ℝ)..L, 2 * Real.sqrt (G' l - G l ^ 2)
      ≤ (2 * Real.sqrt (P * M) / (P + M)) * Real.log ((1 + L / M) / (1 - L / P)) := by
  have := core_strong L P M Lam hL hM hP hstar1 hstar2 G G' hGc hGd hG'c h1 h2 h3
  have hP0 : 0 < P := by linarith
  have : 0 ≤ Real.sqrt (P * M) * G 0 := mul_nonneg (Real.sqrt_nonneg _) h4
  linarith

/-- Strictness for `G 0 > 0` (every smooth generator has `G 0 > 0`). -/
theorem core_strict
    (L P M Lam : ℝ) (hL : 0 < L) (hM : 0 < M) (hP : L < P)
    (hstar1 : L / ((P - L) * (M + L)) = 1)
    (hstar2 : ∫ l in (0:ℝ)..L, l / ((P - l) * (M + l)) = Lam)
    (G G' : ℝ → ℝ)
    (hGc : ContinuousOn G (Icc 0 L))
    (hGd : ∀ x ∈ Ioo 0 L, HasDerivAt G (G' x) x)
    (hG'c : ContinuousOn G' (Icc 0 L))
    (h1 : G L = 1) (h2 : ∫ l in (0:ℝ)..L, G l = Lam)
    (h3 : ∀ l ∈ Icc 0 L, G l ^ 2 ≤ G' l) (h4 : 0 < G 0) :
    ∫ l in (0:ℝ)..L, 2 * Real.sqrt (G' l - G l ^ 2)
      < (2 * Real.sqrt (P * M) / (P + M)) * Real.log ((1 + L / M) / (1 - L / P)) := by
  have := core_strong L P M Lam hL hM hP hstar1 hstar2 G G' hGc hGd hG'c h1 h2 h3
  have hP0 : 0 < P := by linarith
  have : 0 < Real.sqrt (P * M) * G 0 :=
    mul_pos (Real.sqrt_pos.2 (mul_pos hP0 hM)) h4
  linarith

end

end NullCone
