# Reproduce

Python 3 with mpmath, numpy, scipy. One core, from this directory, each under a few minutes:

`python3 -I shear_clock_check.py` (Section 2 identities, Theorem 1 quadratures; log `shear_clock_check.log`)
`python3 -I shear_clock_step2.py` (Theorem 1 closed form; log `shear_clock_step2.log`)
`python3 -I shear_clock_step3.py` (convexity of F, inverse; log `shear_clock_step3.log`)
`python3 -I shear_clock_step4.py` (elementary bound, Table 1, small excess, existence curve, whole cut; log `shear_clock_step4.log`)
`python3 -I shear_clock_deeper.py` (Section 6: concavity of Phi, large-clock law, patch-bound ratios, the family of Section 6.2; log `shear_clock_deeper.log`)
`python3 -I shear_clock_figure.py ../figures/fig1_dmax.pdf` (Figure 1: D_max, elementary and (D.52a) bounds over tau; needs matplotlib)
