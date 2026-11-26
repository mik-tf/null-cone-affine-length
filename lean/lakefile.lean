import Lake
open Lake DSL

package nullcone where
  leanOptions := #[⟨`autoImplicit, false⟩]

-- Mathlib revision (tag v4.36.0).
require mathlib from git "https://github.com/leanprover-community/mathlib4" @ "v4.36.0"

@[default_target] lean_lib NullCone
lean_lib Axioms
