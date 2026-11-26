# Affine length bounds the conformal excursion on a null cone: a sharp focusing inequality

Author: Michaël Perreault

Version 1.0.

On a vacuum null cone, a generator that stays expanding within a given normalised affine length moves its conformal 2-metric, a curve in the hyperbolic plane, by less than an explicit sharp constant. The paper gives the closed form, the unique extremal, a sharp whole-cut inequality, and consequences (a sharp patch bound, a gap in a one-loop check, higher dimensions, matter).

Source: [paper.md](paper.md). Text: CC BY 4.0 ([LICENSE-paper.md](LICENSE-paper.md)). Code in `code/`: MIT ([LICENSE](LICENSE)).

Archive contents: `paper.md` and the PDF; `code/` (six Python scripts and logs, `code/REPRODUCE.md`); `lean/` (Lean 4 with Mathlib: `lean/NullCone/Core.lean`, `lean/README.md`). `lean/` verifies only the core inequality of Theorem 1 (`core_lemma`, `core_strong`, `core_strict`); `lean/wip/` holds unfinished drafts, not built and not part of the verified claim. The `.lake/` build directory is excluded.

## Re-run

Python 3 with mpmath, numpy and scipy; see [code/REPRODUCE.md](code/REPRODUCE.md). Lean: see [lean/README.md](lean/README.md).

## Cite as

Perreault, M. Affine length bounds the conformal excursion on a null cone: a sharp focusing inequality, version 1.0. https://github.com/mik-tf/null-cone-affine-length (DOI to be added). See [CITATION.cff](CITATION.cff).
