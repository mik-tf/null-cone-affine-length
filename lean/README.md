# Lean 4 formalisation of the core inequality (Theorem 1)

Affine length bounds the conformal excursion on a null cone: a sharp focusing inequality.

Lean `v4.36.0`, Mathlib tag `v4.36.0` (pinned in `lake-manifest.json`). Build: `lake exe cache get`, then `lake build`.
`Axioms.lean` prints `#print axioms` for the three theorems below.

## What is proved (no `sorry`, axioms only `propext`, `Classical.choice`, `Quot.sound`)

File `NullCone/Core.lean`. Notation: `L = lam0`, `P = r_+`, `M = r_-`, `Lam = Lambda`,
`Q(l) = (P - l)(M + l)`, `G*(l) = l / Q(l)`, `c = 2 sqrt(P M)`.

| Lean | Paper |
|---|---|
| `core_strong` | the optimality step of the proof of Theorem 1: `ell(G) <= ell(G*) - (c/2) G(0)`, with `ell(G*) = D_max = (2 sqrt(P M)/(P+M)) ln((1+beta)/(1-alpha))`, `alpha = L/P`, `beta = L/M` |
| `core_lemma` | Theorem 1, bound `ell(G) <= D_max` for every admissible `G` (`G(0) >= 0`) |
| `core_strict` | Theorem 1, strict inequality `ell(G) < D_max` when `G(0) > 0` (every smooth generator) |

Admissible `G` means: continuous on `[0, L]`, differentiable on `(0, L)` with continuous derivative `G'` on `[0, L]`
(this is the paper's `C^1[0, L]`, stated without one-sided derivatives), `G(L) = 1`, `int_0^L G = Lam`,
`G' >= G^2`. The pair `(P, M)` enters as the paper's two equations, as hypotheses: `G*(L) = 1`
(`hstar1`, first equation of Theorem 1) and `int_0^L G* = Lam` (`hstar2`, second equation, in the form `int l/Q`).
Also proved on the way: the closed form `int_0^L 2m/Q = (2m/(P+M)) ln((1+L/M)/(1-L/P))`
(`integral_inv_Q`, the paper's partial fractions for `ell(G*)`).

Proof: pointwise tangent bound `c u s <= u^2 q - 2 l u g + l^2 + m^2` (AM-GM plus `(u g - l)^2 >= 0`, lemma `ptwise`),
linear in `(G, G')`; integration by parts of `Q (G - G*)` (the paper's `Q' + 2 l = r_+ - r_-` constant);
`G(L) = G*(L)`, `int G = int G*`.

## What stays on paper

The GR-to-1D reduction (Section 2: Raychaudhuri, the hyperbolic metric, the clock identity, `G = e^(-t-X)`), the
numerics, Lemma 2, Proposition 2, Theorem 3 and Sections 5 to 6.

## Work in progress (not built, not part of the claim)

`wip/Existence.lean` (existence of `(alpha, beta)` by the intermediate value theorem) and `wip/Theorem1.lean` (the core
lemma restated in `(alpha, beta)`, contains `sorry`) are drafts. They do not compile yet and are outside the lake
libraries. Next step is in `../STATE.md`.
