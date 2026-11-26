---
title: "Affine length bounds the conformal excursion on a null cone: a sharp focusing inequality"
author: "Michaël Perreault"
date: "2026-11-26"
format: "article"
toc: false
header_footer_policy: "all"
footer: "© 2026 Michaël Perreault · CC BY 4.0"
pageof: true
date_footer: "YYYY-MM-DD"
header_title: "Affine length bounds the conformal excursion"
---

*Corresponding author: Michaël Perreault, independent researcher, ORCID 0009-0003-0017-7658.*

**Abstract.** A bundle of light rays that keeps expanding can twist the shape of its cross-sections only a bounded amount for a given affine length, and we give the exact bound. Precisely, on a vacuum null cone, let a generator stay expanding from $\varepsilon\rho_S$ to $\rho_S$ (the area density is $\rho^2$), within normalised affine length $\tau\rho_S$. Then its conformal 2-metric, a curve in the hyperbolic plane, moves a hyperbolic distance less than an explicit $D_{\max}(\tau,\varepsilon)$, and this constant is sharp. The affine length is $\rho_S\int e^{\Xi}dv$, where $\Xi$ is a quarter of the curve's energy, as in Chandrasekaran ([arXiv:2611.34924](https://arxiv.org/abs/2611.34924)), whose Appendix D has a non-sharp bound of the same kind. New here, all proved: the closed form of the sharp constant and its unique extremal (speed inversely proportional to the remaining affine length), the small-excess law, the whole-cut inequality $\bar u_*\ge(A/4\pi)^{1/2}F_\varepsilon(\langle d\rangle_A)$, a proof that the logarithmic growth of Chandrasekaran's patch bound is necessary, the extension to dimension $d+2$ (factor $(d/2)^{1/2}$) and to matter, and a gap, with an explicit family, in the one-loop check of its Appendix D.3, not in its finiteness proof. Computed only, and labelled so in the text: the matching large-clock law, which suggests that the prefactor $1/\varepsilon$ of that patch bound is loose by about $9$, $54$ and $370$ at $\varepsilon=0.1$, $0.01$, $0.001$ (for $A^*/A_P$ between $10^3$ and $10^{12}$). At $\varepsilon=10^{-3}$ the sharp bound improves the elementary Raychaudhuri bound by $27$ percent at $\tau=2$ and by a factor $4$ at $\tau=100$ (computed from the closed form).

**Keywords.** null cone; Raychaudhuri equation; focusing; hyperbolic plane; affine parameter; characteristic initial data.

**MSC 2020.** 83C05, 83C30, 53C22, 49K15.

# Introduction {#p8.intro}

Free data for vacuum general relativity on a null cone consist of a conformal 2-metric $\bar q_{AB}$ (unimodular, positive) on each generator, as a function of an area parameter, together with an area density on one cut [Sachs62, Rei12, Rei07]. At each generator the matrix $\bar q$ lies in $SL(2,\mathbb R)/SO(2)$, the hyperbolic plane $\mathbb H^2$ of curvature $-1$. So the conformal metric along a generator is a curve $\gamma$ in $\mathbb H^2$, and its net change between two cuts is a hyperbolic distance $d$, which we call the *conformal excursion*. Its linearisation at null infinity is the linear shear memory [BT87]; we do not call $d$ memory, since the nonlinear memory of [Chr91] is a different effect.

Contrapositively, a generator whose conformal metric moves farther than $D_{\max}$ focuses before affine length $\tau\rho_S$, as in the trapped-surface literature [Chr08, KR09, AL14].

**Prior art and what is new.** The setting and the identity $u_*=\rho_S\int_\varepsilon^1e^{\Xi}dv$ are in Chandrasekaran [Cha26, (3.26) to (3.27)], with $\Xi$ as in (3.26c) and the metric in (4.21b to c). Its Appendix D already bounds a total variation of the conformal metric by the clock: it defines $D(y)=\int_\varepsilon^1[(\partial_v\bar q)_{AB}\partial_v\bar q^{AB}]^{1/2}dv$ in (D.51), equal to $\sqrt2$ times our hyperbolic length $\ell$, and proves by Cauchy-Schwarz and Hölder $D(y)^2\le 8[1-e^{-t\Xi_y(1)}](\varepsilon t)^{-1}\int_\varepsilon^1e^{t\Xi_y}dv$ in (D.52a) for $t\in(0,1]$, with an averaged form (D.52b); (D.53) to (D.54) use convexity of the clock along geodesic interpolation, which overlaps our Lemma 2. They serve another purpose there (small symplectic area of local deformations) and are not sharp.

*New here:* (i) the exact best constant $D_{\max}(\tau,\varepsilon)$ in closed form, with its unique extremal and a proof that it is a supremum not attained (Theorem 1); (ii) the small-excess law (Proposition 2); (iii) the sharp whole-cut inequality with convex $F_\varepsilon$ (Theorem 3), which contains the classical area bound as its case $d=0$ and improves it whenever the excursion is nonzero. *Consequences (Section 6):* a sharp patch bound, for which the logarithmic growth of (D.52b) is proved necessary and the looseness of its prefactor $1/\varepsilon$ is computed (by a factor of about $e/(\varepsilon\ln(1/\varepsilon))$ up to a slowly varying factor, over the range tested); the observation that the clock alone cannot bound the shear energy used in D.3, a gap in a one-loop consistency check and not in the finiteness proof of Section 5.2 of [Cha26]; the extension to every dimension (factor $(d/2)^{1/2}$) and to matter satisfying the null energy condition. *Not new:* the clock identity, the use of convexity along geodesics, and the existence of some length-clock bound.

Section 2 re-derives [Cha26, (3.26)] in classical terms. Section 3 proves Theorem 1, Section 4 the whole-cut theorem, Section 5 compares with other bounds, Section 6 gives the consequences, Section 7 the checks and related work.

# Setup and the clock identity {#p8.setup}

*Hypotheses, used throughout.* Vacuum; smooth characteristic data on a cone; every generator expanding from cap to cut, $\theta>0$ (so the area can serve as parameter); the cap is defined generator by generator by $\rho=\varepsilon\rho_S(y)$ as in [Cha26, (3.18)]; the affine parameter $u$ is normalised by $p_\varepsilon:=\partial_u\rho|_{\rm cap}=1$ (the clock of [Cha26, (3.27)] is $p_\varepsilon\Delta u$, which is invariant under rescaling of $u$ and equals $\Delta u$ in this normalisation); the reference area form $\epsilon_0$ has total $\int\epsilon_0=4\pi$. The cap is a regulator inherited from [Cha26]. 

Fix a generator and an affine parameter $u$. Write $q_{AB}$ for the transverse metric and $\rho=(\det q)^{1/4}$ (relative to $\epsilon_0$), so that the area element is $\rho^2\epsilon_0$ and the area density is $\rho^2$; put $\bar q=q/\rho^2$, which has unit determinant, as in [Cha26, (3.17c)]. Since $q^{-1}q'=2(\rho'/\rho)\mathbb 1+\bar q^{-1}\bar q'$ with trace-free second term (prime is $d/du$), the expansion is $\theta=\tfrac12\operatorname{tr}(q^{-1}q')=2\rho'/\rho$, the shear is $\tfrac12\bar q^{-1}\bar q'$ (as a mixed tensor), $\sigma_{AB}\sigma^{AB}=\tfrac14\operatorname{tr}((\bar q^{-1}\bar q')^2)$, and the Raychaudhuri equation $\theta'=-\tfrac12\theta^2-\sigma_{AB}\sigma^{AB}$ reads $2\rho''/\rho-2\rho'^2/\rho^2=-2\rho'^2/\rho^2-\sigma_{AB}\sigma^{AB}$, that is
$$
\rho''=-\tfrac12\,\sigma_{AB}\sigma^{AB}\,\rho .
\tag{2.1}
$$

**The metric.** On unimodular positive symmetric matrices the form $g=\tfrac12\operatorname{tr}(\bar q^{-1}d\bar q\,\bar q^{-1}d\bar q)$ is the hyperbolic metric of curvature $-1$: for $\bar q=\operatorname{diag}(e^{z_1},e^{-z_1})$ it is $dz_1^2$, and on $ac-b^2=1$ it is $-dx_0^2+dx_1^2+dx_2^2$ with $x_0=(a+c)/2$, $x_1=(a-c)/2$, $x_2=b$ (compare [Cha26, (4.21b to c)]). The distance satisfies
$$
\cosh d(\bar q_0,\bar q_1)=\tfrac12\operatorname{tr}(\bar q_0^{-1}\bar q_1).
\tag{2.2}
$$
With $|\gamma_u|_g$ the speed, $\sigma_{AB}\sigma^{AB}=\tfrac12|\gamma_u|_g^2$ and (2.1) reads $\rho''=-\tfrac14|\gamma_u|_g^2\rho$.

**The clock.** Take $v=\rho/\rho_S$ as parameter: the cut is $v=1$, the cap $v=\varepsilon$. With $p=d\rho/du$, $|\gamma_u|^2=p^2|\gamma_\rho|^2$ and (2.1) gives $dp/d\rho=-\tfrac14p\rho|\gamma_\rho|^2$, so $p(v)=p_\varepsilon e^{-\Xi(v)}$ with
$$
\Xi(v)=\tfrac14\int_\varepsilon^v s\,|\partial_s\gamma|_g^2\,ds
\tag{2.3}
$$
([Cha26, (3.26a to c), (4.21b to c)]; at the cut $\Xi(1)=\log(p_\varepsilon/\partial_u\rho|_S)$, [Cha26, (C.16b)]). Since $du=\rho_S\,dv/p$,
$$
u_*(y)=\rho_S(y)\int_\varepsilon^1 e^{\Xi(v)}\,dv .
\tag{2.4}
$$
With $w=\ln s$ one has $s|\gamma_s|^2ds=|\gamma_w|^2dw$, so $\Xi(v)=\tfrac14\int_{\ln\varepsilon}^{\ln v}|\gamma_w|_g^2dw=\tfrac14E[\gamma]$, a quarter of the energy of $\gamma$ in $\mathbb H^2$ (this follows from (2.3) and the substitution, not from (2.2)). Put $t=\ln(v/\varepsilon)\in[0,\Lambda]$, $\Lambda=\ln(1/\varepsilon)$, $a(t)=|\gamma_t|_g\ge0$ and $X(t)=\Xi$. Then $X'=a^2/4$, $X(0)=0$ and
$$
\tau:=\frac{u_*}{\rho_S}=\varepsilon\int_0^\Lambda e^{t+X(t)}\,dt .
\tag{2.5}
$$
The clock depends on the curve only through its speed $a(t)$; the hyperbolic length $\ell=\int_0^\Lambda a\,dt$ bounds the excursion, $d=d(\gamma(0),\gamma(\Lambda))\le\ell$. Always $\tau\ge1-\varepsilon$, with equality iff $a\equiv0$.

# The per-generator theorem {#p8.gen}

For $D\ge0$ let $F_\varepsilon(D)=\inf\{\tau(a):\ a\ge0,\ a\in L^2,\ \int_0^\Lambda a\,dt=D\}$ ($F_\varepsilon(0)=1-\varepsilon$; the infimum is not attained among smooth profiles) and $D_{\max}(\tau,\varepsilon)=\sup\{\int a\,dt:\ \tau(a)=\tau\}$, the inverse of $F_\varepsilon$ (Lemma 2).

**Theorem 1.** Under the hypotheses of Section 2, let $\tau>1-\varepsilon$. There are unique $\alpha\in(0,1)$, $\beta>0$ with
$$
\frac\tau\varepsilon=\frac{\alpha\beta}{(1-\alpha)(1+\beta)},\qquad
\ln\frac1\varepsilon=\frac{\beta\ln\frac1{1-\alpha}-\alpha\ln(1+\beta)}{\alpha+\beta},
$$
and with them
$$
D_{\max}(\tau,\varepsilon)=\frac{2\sqrt{\alpha\beta}}{\alpha+\beta}\,\ln\frac{1+\beta}{1-\alpha}.
$$
Every smooth generator with clock $\tau$ has hyperbolic length $\ell<D_{\max}(\tau,\varepsilon)$, hence excursion $d<D_{\max}$. The constant is a supremum over smooth data, not attained, approached only along sequences converging to the profile $a=c/\lambda$, where $\lambda(t)=\int_t^\Lambda e^{s+X(s)}ds$ is the remaining clock and $c$ a constant (speed inversely proportional to the remaining clock). In parametric form it is $t(\lambda)=\Lambda-\frac1{r_++r_-}\bigl[r_+\ln\frac{r_+}{r_+-\lambda}-r_-\ln\frac{r_-+\lambda}{r_-}\bigr]$ ($r_\pm$ below); it has $\Xi(1)=\infty$ and is not admissible.

*Proof.* The remaining clock decreases strictly from $\lambda_0=\tau/\varepsilon$ to $0$. Use $\lambda$ as variable and $G(\lambda)=e^{-t-X}$, smooth on $[0,\lambda_0]$ because $d\lambda/dt=-e^{t+X}\ne0$. The chain rule gives $G'=(1+X')G^2$, and
$$
\int_0^{\lambda_0}G\,d\lambda=\Lambda,\qquad G(\lambda_0)=1,\qquad G'\ge G^2,\qquad
\ell=2\int_0^{\lambda_0}\sqrt{G'-G^2}\,d\lambda ,
$$
since $a=2\sqrt{X'}$, $X'=G'/G^2-1$, $dt=-G\,d\lambda$. Smooth data have $G(0)=\varepsilon e^{-\Xi(1)}>0$. Let $\mathcal K$ be the set of $G\in C^1[0,\lambda_0]$ satisfying these constraints with $G(0)\ge0$. It is convex, since $G\mapsto G^2$ is convex. The function $h(g,q)=\sqrt{q-g^2}$ is strictly concave on $\{q\ge g^2\}$ ($\sqrt\cdot$ is concave increasing, $q-g^2$ concave, and a segment on which $h$ is affine forces $dg=dq=0$), so $\ell$ has at most one maximiser on $\mathcal K$.

*Candidate.* Let $Q(\lambda)=(r_+-\lambda)(r_-+\lambda)=\tfrac{c^2}4+(r_+-r_-)\lambda-\lambda^2$, $c^2=4r_+r_-$, and $G^*=\lambda/Q$. Then $G^{*\prime}-G^{*2}=c^2/(4Q^2)$, so $a=c/\lambda$ (as $X'=c^2/4\lambda^2$). Put $\alpha=\lambda_0/r_+$, $\beta=\lambda_0/r_-$. Then $G^*(\lambda_0)=1$ is the first equation of the theorem; partial fractions $\lambda/Q=\frac1{r_++r_-}[\frac{r_+}{r_+-\lambda}-\frac{r_-}{r_-+\lambda}]$ turn $\int G^*=\Lambda$ into the second; and $\ell(G^*)=c\int_0^{\lambda_0}d\lambda/Q=\frac c{r_++r_-}\ln\frac{1+\beta}{1-\alpha}$ is the stated $D_{\max}$ (indeed $c/(r_++r_-)=2\sqrt{\alpha\beta}/(\alpha+\beta)$).

*Existence of $(\alpha,\beta)$.* The first equation gives $\beta=r/(1-r)$, $r=\lambda_0(1-\alpha)/\alpha$, for $\alpha\in(\lambda_0/(1+\lambda_0),1)$. Along this curve the right side $\Lambda(\alpha)$ of the second equation is continuous, tends to $\ln(1+\lambda_0)$ as $\alpha\downarrow\lambda_0/(1+\lambda_0)$ ($\beta\to\infty$) and to $0$ as $\alpha\to1$ ($\beta\to0$). Since $\tau>1-\varepsilon$ means $1+\lambda_0>1/\varepsilon$, i.e. $\ln(1+\lambda_0)>\ln(1/\varepsilon)$, the intermediate value theorem gives a solution. For it $r_+>\lambda_0$, so $Q>0$ on $[0,\lambda_0]$ and $G^*\in\mathcal K$, with $G^*(0)=0$.

*Optimality and strictness.* On $[0,\lambda_0]$, $G^{*\prime}-G^{*2}\ge c^2/(4\max Q^2)>0$, so $h$ is smooth near $(G^*,G^{*\prime})$. For $G\in\mathcal K$ put $\eta=G-G^*$; $G_s=G^*+s\eta\in\mathcal K$ and $s\mapsto\ell(G_s)$ is concave, so $\ell(G)\le\ell(G^*)+\frac d{ds}\ell(G_s)|_{0^+}$. Using $h_g=-2\lambda/c$ and $h_q=Q/c$ at $G^*$ and integrating by parts with $\eta(\lambda_0)=0$, $\int\eta=0$,
$$
\frac{d}{ds}\ell(G_s)\Big|_{0^+}=\frac2c\Bigl[Q\eta\Bigr]_0^{\lambda_0}-\frac2c\int_0^{\lambda_0}\bigl(Q'+2\lambda\bigr)\eta\,d\lambda=-\frac c2\,G(0)\le0,
$$
because $Q'+2\lambda=r_+-r_-$ is constant. Hence $\ell(G)\le\ell(G^*)-\frac c2G(0)$, with strict inequality for smooth data ($G(0)>0$). So $G^*$ is the unique maximiser on $\mathcal K$. This also gives uniqueness of $(\alpha,\beta)$: two solutions give two maximisers, and $G^*=\lambda/Q$ determines $Q$, hence $(r_+,r_-)$.

*Sharpness.* Take a smooth feasible $G_0$ with the same $(\tau,\varepsilon)$ and $G_s=(1-s)G^*+sG_0$: feasible, smooth, $G_s(0)>0$, $G_s'-G_s^2\ge(1-s)c^2/(4Q^2)>0$, and $\ell(G_s)\ge(1-s)\ell(G^*)+s\ell(G_0)\to\ell(G^*)$ as $s\to0$. $\square$

Near its end $e^\Xi\sim c^2/(4\lambda)$ and $\lambda\sim(\Lambda-t)^{1/2}$, so $\Xi\to\infty$ logarithmically while the clock converges; $\theta_\ell=0$ at the cut, so a bounded clock does not bound $\Xi(1)$.

**Lemma 2.** $F_\varepsilon$ is convex, continuous at $0$ and strictly increasing on $[0,\infty)$; hence $D_{\max}=F_\varepsilon^{-1}$ is well defined on $(1-\varepsilon,\infty)$.

*Proof.* The length is affine in $a$, $X_a=\tfrac14\int a^2$ is convex in $a$ and $e^{t+X}$ is convex increasing in $X$, so $\tau(a)$ is convex and its infimum under an affine constraint is convex in $D$. With $F_\varepsilon(0)=1-\varepsilon<F_\varepsilon(D)$ for $D>0$ (Theorem 1) it is strictly increasing; continuity at $0$ is Proposition 2. $\square$

**Proposition 2 (small excess).** Let $\delta=F_\varepsilon(D)-(1-\varepsilon)$. As $D\to0$ at fixed $\varepsilon$, $\beta\to\infty$ and, with $B=\ln(1/\varepsilon)+\ln\beta$,
$$
D^2=\frac{4(1-\varepsilon)B^2}{\beta}(1+o(1)),\qquad \delta=\frac{(1-\varepsilon)(B-1)}{\beta}(1+o(1)),\qquad\text{so}\qquad \delta=\frac{D^2(B-1)}{4B^2}(1+o(1)),
$$
and $B=\ln\frac1{\varepsilon D^2}+2\ln B+\ln(4(1-\varepsilon))+o(1)$. In particular $\delta\sim D^2/(4\ln(1/(\varepsilon D^2)))$ (the relative error is $O(\ln\ln(1/D)/\ln(1/D))$), and $D^2\approx4\delta\bigl(\ln\frac1\varepsilon+\ln\frac1\delta\bigr)$ to leading order, with relative error of order $1/\ln(1/\delta)$ (computed ratio $1.1$ to $1.3$ for $\delta$ from $10^{-4}$ to $10^{-14}$).

*Proof.* A finite limit point of $\beta$ along $D\to0$ would give, at $\tau=1-\varepsilon$, a solution of the two equations with $D_{\max}>0$, whereas clock $1-\varepsilon$ forces $a\equiv0$; and $\beta\to0$ is excluded because $\Lambda(\alpha)\to0\ne\ln(1/\varepsilon)$. So $\beta\to\infty$. The first equation gives $1-\alpha=\beta/((1+\beta)\lambda_0+\beta)$, so with $1+\lambda_0=(1+\delta)/\varepsilon$ and $\kappa=\lambda_0/(1+\lambda_0)\to1-\varepsilon$, $\ln\frac1{1-\alpha}=\Lambda+\ln(1+\delta)+\kappa/\beta+O(\beta^{-2})$. The second equation, multiplied by $\alpha+\beta$, becomes $\beta\ln(1+\delta)=\alpha(\Lambda+\ln(1+\beta))-\kappa+o(1)$, and $\alpha\to\kappa$, giving $\beta\delta=\kappa(B-1)+o(1)$. In $D$, $\ln\frac{1+\beta}{1-\alpha}=B+o(1)$ and $2\sqrt{\alpha\beta}/(\alpha+\beta)=2\sqrt{\kappa/\beta}(1+o(1))$, so $D^2=4\kappa B^2/\beta\,(1+o(1))$. Dividing gives $\delta/D^2=(B-1)/(4B^2)(1+o(1))$; solving $\beta=4\kappa B^2/D^2$ for $\ln\beta$ gives the implicit equation for $B$. $\square$

Numerically $\delta$ matches the parametric form to $1.0000$ ($\varepsilon=10^{-1}$ to $10^{-12}$, $\beta$ up to $10^{10}$); the leading form $D^2/(4\ln(1/(\varepsilon D^2)))$ is approached only logarithmically (ratio $0.42$ to $0.80$ in that range), so the implicit form is the useful statement.

# The whole-cut theorem {#p8.cut}

On a cut $S$ with area density $\rho_S$ and area $A=\int_S\rho_S^2\epsilon_0$ let $d(y)=d(\bar q(\varepsilon,y),\bar q(1,y))$, by (2.2). Let $\langle d\rangle_A=\int\rho_S^2d\,\epsilon_0/A$ and $\bar u_*=\int u_*\rho_S\epsilon_0/\int\rho_S\epsilon_0$.

**Theorem 3.** Under the hypotheses of Section 2, with a smooth cut,
$$
\bar u_*\ \ge\ \Bigl(\frac A{4\pi}\Bigr)^{1/2}F_\varepsilon\bigl(\langle d\rangle_A\bigr),
$$
and the same with $d(y)$ replaced by the hyperbolic length of $\gamma$. The inequality is sharp: equality holds in the supremum sense when $\rho_S$ is constant and every generator approaches the extremal of Theorem 1 of the same length.

*Proof.* By Theorem 1, $u_*(y)\ge\rho_S(y)F_\varepsilon(d(y))$ (as $F_\varepsilon$ is increasing and $d\le\ell$). Integrating against $\rho_S\epsilon_0$ and using Jensen and Lemma 2, $\bar u_*\int\rho_S\epsilon_0\ge A\langle F_\varepsilon(d)\rangle_A\ge A\,F_\varepsilon(\langle d\rangle_A)$; and $\int\rho_S\epsilon_0\le(4\pi A)^{1/2}$ by Cauchy-Schwarz. $\square$

The per-generator problems decouple (no angular derivatives in the clock; free conformal data).

**Classical corollary.** Since $F_\varepsilon(0)=1-\varepsilon$ and $F_\varepsilon$ is increasing, $\bar u_*\ge(1-\varepsilon)(A/4\pi)^{1/2}$ for all data: the classical statement behind footnote 8 and (5.21c) of [Cha26] (an operator inequality there, which we do not recover). Every classical solution with $\bar u_*\le L$ satisfies $A\,F_\varepsilon(\langle d\rangle_A)^2\le4\pi L^2$; at $\langle d\rangle_A=0$ this is the classical content of (5.21d).

# Comparison with other bounds {#p8.comp}

*Elementary bound.* Integrating (2.1) from cap to cut, $p_S-p_\varepsilon=-\tfrac14\int|\gamma_u|^2\rho\,du$, so $\int|\gamma_u|^2\rho\,du\le4p_\varepsilon=4$. By Cauchy-Schwarz, $\ell^2\le4\int du/\rho$. Since $\rho$ is concave in $u$ it lies above its chord from $\varepsilon\rho_S$ to $\rho_S$ over $\Delta u=\tau\rho_S$, so $\int du/\rho\le\tau\ln(1/\varepsilon)/(1-\varepsilon)$, and
$$
\ell\le L_{\rm el}(\tau,\varepsilon)=2\Bigl(\frac{\tau\ln(1/\varepsilon)}{1-\varepsilon}\Bigr)^{1/2}.
$$
It grows like $\tau^{1/2}$ and does not vanish as $\tau\to1-\varepsilon$.

*Appendix D of [Cha26].* With $\int_\varepsilon^1e^{t\Xi}dv\le\tau^t(1-\varepsilon)^{1-t}$ (Hölder), $D=\sqrt2\ell$, and the factor $1-e^{-t\Xi_y(1)}$ bounded by $1$ (the only bound the clock gives, since it does not bound $\Xi(1)$), (D.52a) gives $\ell^2\le(4/\varepsilon)\min_{t\in(0,1]}\tau^t(1-\varepsilon)^{1-t}/t$. The factor $1/\varepsilon$ makes it weak here; it is designed for small patches.

**Table 1** ($\varepsilon=10^{-3}$, exact closed form, 30 digits). Largest hyperbolic length compatible with clock $\tau$:

| $\tau$ | $D_{\max}$ (sharp) | $L_{\rm el}$ | $L_{\rm el}/D_{\max}$  | (D.52a), bound on $\ell$, not $D$ |
|---|---|---|---|---|
| 1.01 | 0.812 | 5.29 | 6.5 | 63.6 |
| 2 | 5.414 | 7.44 | 1.37 | 89.4 |
| 10 | 9.269 | 16.6 | 1.79 | 158 |
| 100 | 12.65 | 52.6 | 4.16 | 224 |

Figure 1 shows the three bounds over $\tau$. The sharp bound is below the elementary one for every $\tau$ (as it must be) and the gain is smallest near $\tau=2$ (the elementary bound is $37$ percent above $D_{\max}$ there, $D_{\max}$ is $27$ percent below it); the sharp bound vanishes as $\tau\to1-\varepsilon$ and for large $\tau$ it replaces growth $\tau^{1/2}$ by growth like $\sqrt{4\ln(1/\varepsilon)\ln\tau}$ (lower bound proved in Section 6.1, the matching upper bound computed there). At fixed $\tau$, $D_{\max}\to\infty$ as $\varepsilon\to0$ ($3.75$, $4.71$, $6.17$, $8.30$ at $\tau=1.5$, $\varepsilon=10^{-2},10^{-4},10^{-8},10^{-16}$), so the bound needs the cap. Whether vertex regularity (bounded $|\partial_s\gamma|$ as $s\to0$) gives an $\varepsilon$-uniform bound is open.

![The sharp bound $D_{\max}$, the elementary bound $L_{\rm el}$ and the Appendix D bound (D.52a) over the clock excess $\tau-(1-\varepsilon)$, $\varepsilon=10^{-3}$ (computed from the closed form).](figures/fig1_dmax.pdf){width=75%}

At fixed endpoints the clock minimiser lies on the geodesic: rescaling the speed of a given curve by $d/\ell$ (so that it runs along the geodesic with the same parametrisation) lowers $\Xi$ and keeps the endpoint.

# Consequences and extensions {#p8.ext}

This section collects four consequences of Theorem 1. Each is labelled *proved* or *computed*, and where a step rests on a numerical check this is said at the step.

## A sharp patch bound {#p8.patch}

In the one-loop analysis of [Cha26, Appendix D.3] a clock budget on a patch $P$ of area $A_P$ is combined with the total-variation bound (D.52a) to bound $\int_PD^2\,dA$, where $D=\sqrt2\,\ell$ (D.51). Write the budget as in (D.50): $\int_P\tau_y\,dA\le(1-\varepsilon)A^*$. Let $\Phi(\tau)=D_{\max}(\tau,\varepsilon)^2$ for $\tau>1-\varepsilon$ and $\Phi(1-\varepsilon)=0$. It is increasing and $\Phi\le L_{\rm el}^2$, which is linear in $\tau$. Let $\hat\Phi$ be its least concave majorant on $[1-\varepsilon,\infty)$; it is finite, increasing and concave.

**Proposition 4 (proved).** Under the budget (D.50),
$$
\int_PD(y)^2\,dA\ \le\ 2A_P\,\hat\Phi\!\Bigl(\frac{(1-\varepsilon)A^*}{A_P}\Bigr).
$$
It is sharp at every clock value $x$ with $\hat\Phi(x)=\Phi(x)$: take every generator in $P$ with clock $x=(1-\varepsilon)A^*/A_P$ and close to the extremal of Theorem 1.

*Proof.* By Theorem 1, $D(y)^2=2\ell_y^2\le2\Phi(\tau_y)\le2\hat\Phi(\tau_y)$. Jensen's inequality for the concave function $\hat\Phi$ and the probability measure $dA/A_P$ gives $\int_P\hat\Phi(\tau_y)dA\le A_P\hat\Phi(\langle\tau\rangle_P)$, and $\hat\Phi$ is increasing with $\langle\tau\rangle_P\le(1-\varepsilon)A^*/A_P$. For sharpness, a constant clock $x$ with extremal-like generators gives $2A_P\Phi(x)$ in the supremum sense. $\square$

The proof does not use concavity of $\Phi$ itself. Numerically $\Phi$ is concave: all second divided differences are negative on a grid of 49 points in $\ln\tau$, for $\varepsilon=10^{-1},10^{-3},10^{-6}$ (computed), so numerically $\hat\Phi=\Phi$ there and the bound reads $2A_P\,D_{\max}((1-\varepsilon)A^*/A_P,\varepsilon)^2$. An analytic proof of concavity is not given here.

*Comparison with (D.52b).* The bound (D.52b) of [Cha26] reads $\int_PD^2dA\le8e(1-\varepsilon)(A_P/\varepsilon)[1+\ln(A^*/A_P)]$. Two statements about large clock, with different status.

*Proved (lower half).* For every $\tau>1$, $\Phi(\tau)\ge4\ln(1/\varepsilon)\ln\tau$. Indeed the constant-speed curve $a\equiv2\sqrt{k-1}$, $k>1$, is admissible, has $X=(k-1)t$, clock $\varepsilon(e^{k\Lambda}-1)/k$ and length $2\Lambda\sqrt{k-1}$. Choose $k$ so that the clock equals $\tau$ (possible: the clock increases continuously in $k$, from $1-\varepsilon$ as $k\to1$ to infinity). Then $\tau=\varepsilon(e^{k\Lambda}-1)/k$. From $e^{k\Lambda}>\tau k/\varepsilon$ one gets $k-1>(\ln\tau+\ln k)/\Lambda>\ln\tau/\Lambda$, so the length squared exceeds $4\Lambda\ln\tau$, and $\Phi\ge$ that by Theorem 1. Hence no bound on $\int_PD^2dA$ can grow more slowly than $A_P\ln(A^*/A_P)$ with this coefficient: the logarithm of (D.52b) is necessary, not an artefact.

*Computed, not proved (upper half).* The ratio of $\Phi$ to $4\ln(1/\varepsilon)\ln\tau$ is $1.26$, $1.18$, $1.12$, $1.08$ at $\tau=10^2,10^4,10^8,10^{16}$ ($\varepsilon=10^{-3}$ only), decreasing slowly toward $1$, while the difference $\Phi-4\ln(1/\varepsilon)\ln\tau$ still grows (from $46$ to $77$ over that range). So the large-clock law $\Phi(\tau)\sim4\ln(1/\varepsilon)\ln\tau$ is a numerical observation over $\tau$ up to $10^{16}$ at one $\varepsilon$, and we have no proof of the upper half. Under that law the prefactor of (D.52b) is the loose part: it comes from bounding the weight $1/v$ on $[\varepsilon,1]$ by $1/\varepsilon$, whereas the computed factor is $\ln(1/\varepsilon)$. For $A^*/A_P$ between $10^3$ and $10^{12}$ the ratio of the bound (D.52b) to the sharp bound is about $9$, $54$ and $370$ at $\varepsilon=0.1$, $0.01$ and $0.001$ (computed: ranges $9.1$ to $9.8$, $53$ to $55$, $369$ to $373$; outside that range it differs, for example $417$ at $A^*/A_P=10$, $\varepsilon=10^{-3}$), and it grows without bound as the cap is removed. This is not a criticism of (D.52b), which was written for another purpose.

## The clock does not bound the shear energy of Appendix D.3 {#p8.gap}

In [Cha26, D.3] the one-loop argument bounds $\log d(s)$ after (D.46c) by a quantity $Q^*$, and the sentence following (D.46c) states $Q^*\lesssim_\varepsilon1+[\ell(D)/\rho_S-(1-\varepsilon)]^{1/3}$ for the expanding cone in flat spacetime. Here $L_D:=\ell(D)$ is the clock cap of [Cha26], not the hyperbolic length $\ell$ of this paper. As we read it, $Q[h]$ of (D.40) contains the term $\tfrac18\int v^2|\partial_vh|^2dv$. The normalisation of $h$ in (D.40) against $\gamma$ is not fixed by the source and only changes constants. We observe that the clock alone does not control this term.

*What the clock controls.* Since $e^{\Xi}-1\ge\Xi$ and $\int_\varepsilon^1\Xi\,dv=\tfrac14\int_\varepsilon^1s(1-s)|\gamma_s|^2ds$ (swap the order of integration), the clock excess satisfies
$$
\tau-(1-\varepsilon)\ \ge\ \tfrac14\int_\varepsilon^1s(1-s)\,|\partial_s\gamma|_g^2\,ds ,
$$
a shear energy with weight $s(1-s)$, which degenerates at the cut $s=1$ (proved). Together with Theorem 1 this is all the clock gives: the hyperbolic length and the $(1-s)$-weighted energy.

*The family.* Fix $K>0$ and put $\eta=e^{-K/4}/K$. Let $|\gamma_v|_g^2=K/\eta$ on $[1-\eta,1]$ and $0$ elsewhere (take $\rho_S=1$). Then $\Xi\le K/4+o(1)$, $\int(e^\Xi-1)dv\approx(4\eta/K)e^{K/4}=4/K^2$, so the clock excess is $4/K^2\,(1+o(1))$, while $\tfrac18\int v^2|\gamma_v|^2dv\approx K/8$ and the hyperbolic length is $\sqrt{K\eta}=e^{-K/8}$. Exact quadrature (computed) gives clock excess $0.049,\,0.028,\,0.0025,\,4\times10^{-4},\,10^{-4}$ at $K=1,10,40,100,200$, with $\tfrac18\int v^2|\gamma_v|^2=0.053,\,1.24,\,5.0,\,12.5,\,25.0$. So the clock excess tends to zero while the shear term grows without bound. This is the phenomenon of Section 3: a bounded clock does not bound $\Xi(1)$.

*What this does and does not say.* This is a gap in a sentence of a one-loop consistency check. It is not a counterexample to the finiteness proof of [Cha26, Section 5.2], which, as we read it, uses only the zero-shear area bound (5.21d) and not (D.51) to (D.54). Nor do we claim that the sentence fails: the bound on $\log d(s)$ can still hold through the separately fixed incoming-expansion datum at the cap, which (D.38) names, and not through $\ell(D)$. At fixed radial resolution the family is excluded (the ratio of the shear term to the clock excess stays bounded by a constant times the number of radial modes), which matches the regulator of [Cha26]. We have not checked how the incoming-expansion datum enters the whole argument, and we may have misread the role of (D.46c). One more caveat: (D.40) is a second-order expression in a small perturbation, and the family has amplitude $e^{-K/8}$ (small, the hyperbolic length above), so it lies inside that regime.

## Higher dimensions {#p8.dim}

In spacetime dimension $d+2$ let $r=(\det q)^{1/(2d)}$ (for $d=2$ this is $\rho$). With $\theta=d\,r'/r$ the Raychaudhuri equation $\theta'=-\theta^2/d-\sigma_{AB}\sigma^{AB}$ becomes $r''=-(\sigma_{AB}\sigma^{AB}/d)\,r$. On $SL(d,\mathbb R)/SO(d)$ with the metric $g=\tfrac12\operatorname{tr}(\bar q^{-1}d\bar q)^2$ one has $\sigma_{AB}\sigma^{AB}=\tfrac12|\gamma_u|_g^2$ exactly as for $d=2$, so $r''=-(|\gamma_u|_g^2/2d)\,r$ and the clock identity holds with $X'=a^2/(2d)$. The clock depends on the curve only through its speed $a$. With $a=(d/2)^{1/2}\tilde a$ the problem is the case $d=2$ in the variable $\tilde a$, so (proved)
$$
\ell\ \le\ (d/2)^{1/2}\,D_{\max}(\tau,\varepsilon),\qquad \varepsilon=\text{ratio of }r\text{ at the cap and cut}.
$$
The excursion is bounded by $\ell$ because $SL(d,\mathbb R)/SO(d)$ is a complete simply connected manifold of nonpositive curvature. The bound is sharp: run the extremal speed of Theorem 1, rescaled, along a unit-speed geodesic, which minimises distance in nonpositive curvature. Higher rank is invisible because only arc length enters. The whole-cut inequality and the patch bound follow the same lines with the factor $d/2$ and with the exponents of the area changed; we do not write them out.

## Matter {#p8.matter}

With a stress tensor, $r''=-\tfrac1d(\sigma_{AB}\sigma^{AB}+8\pi G\,T_{uu})\,r$. In the area parameter $v$ the clock exponent splits as $\Xi=\Xi_\gamma+\Xi_m$, with $\Xi_m(v)=\tfrac{8\pi G}{d}\int_\varepsilon^v sT_{ss}\,ds$ (units of $r_S$, $T_{ss}=T(\partial_s,\partial_s)$). Under the null energy condition $\Xi_m\ge0$ and is nondecreasing. Then $e^{\Xi_\gamma+\Xi_m}\ge e^{\Xi_\gamma}(1+\Xi_m)\ge e^{\Xi_\gamma}+\Xi_m$, and integrating (proved)
$$
\tau\ \ge\ F_\varepsilon(\ell)+W,\qquad W=\int_\varepsilon^1\Xi_m\,dv=\int_\varepsilon^1(1-v)\,d\Xi_m(v),
$$
so $\ell\le D_{\max}(\tau-W,\varepsilon)$: the conformal excursion and the $(1-v)$-weighted null energy share one clock budget. The bound is sharp in each variable separately: $W=0$ is Theorem 1, and for $\ell=0$ it is exact to first order in the matter. It is not jointly sharp, since the step $e^{\Xi_\gamma}\ge1$ is loose at large $\ell$; the joint optimum is not solved here. Matter near the cut costs no clock, by the weight $1-v$, the same degeneracy as for $\Xi(1)$.

We do not treat the quantum problem. A bound of the form $\hat\tau-F_\varepsilon(\hat\ell)\ge-C\hbar$ for all states would need control of the ordering of $e^{\hat\Xi}$ against $\hat\rho_S$, and the extremal has $\Xi(1)=\infty$; we have no result.

# Checks and related work {#p8.checks}

Six scripts in `code/` reproduce Section 2, Theorem 1 (quadratures to $10^{-15}$), convexity of $F_\varepsilon$, Table 1, Proposition 2 and the computed statements of Section 6 and Figure 1 (`code/REPRODUCE.md`). The sharpness statements (Theorems 1 and 3, Proposition 4) are over characteristic data with a free conformal curve per generator, smooth at the cone; we do not address existence of the vacuum development up to the cut, and the extremal has $\Xi(1)=\infty$ (stated above). The operator analogue of Theorem 3 is open (no operator inequality is known to us).

Christodoulou [Chr08], Klainerman and Rodnianski [KR09] and An and Luk [AL14] use $\int|\hat\chi|^2$ against the expansion to prove trapped-surface formation from small conformal data; they state no inequality between affine length and hyperbolic excursion (our extremal with zero expansion at the cut is the threshold case). Wieland [Wie24, Wie25] bounds shear through the $SU(1,1)$ structure, a different quantity; [CFL23] uses the Beltrami parametrisation. Characteristic initial value problems for the Einstein equations are treated in [Rend90, Luk12, CCM11]. The trapped-surface results of [KLR14] extend [Chr08]. The nonlinear memory effect is that of [Chr91]. We found no statement of Theorems 1 or 3 beyond [Cha26, Appendix D].

**Formal verification.** The core inequality of Theorem 1 (the pointwise bound, its integrated form $\ell(G)\le D_{\max}$ and its strict form for $G(0)>0$; in the Lean files `core_lemma`, `core_strong` and `core_strict`) is formally verified in Lean 4 with Mathlib, using only the standard axioms (`propext`, `Classical.choice`, `Quot.sound`). The files are in the directory `lean/` of the paper's code and data archive. The existence of $(\alpha,\beta)$ and the sharpness rest on the written proof; the reduction from the null cone to the one-dimensional problem and the numerics are not formalised. Files in `lean/wip/` are drafts, not built, and not part of the verified claim.

# Code and data {#p8.code .unnumbered}

The six Python scripts with their logs, the reproduction notes and the Lean 4 formalisation of the core inequality of Theorem 1 are archived at <https://github.com/mik-tf/null-cone-affine-length> and on Zenodo, <https://doi.org/10.5281/zenodo.24572220> (concept DOI, all versions). `code/REPRODUCE.md` explains how to reproduce each result.

# Acknowledgements {#p8.ack .unnumbered}

We are grateful to Jason Glowney, whose geometric treatment of light, the light cone and hyperbolic curvature in *Arithmoí*, Book 1 [Glo23], led us to the hyperbolic-plane view of null data on which this paper is built.

# Use of AI {#p8.ai .unnumbered}

Claude (Anthropic) was used in November 2026, under the authors' direction, for exploratory computation, derivations, code, checks, literature search and drafting.

# References {#p8.refs .unnumbered}

An asterisk marks a source not read in full; statements about it rest on abstracts or the passages named (CFL23 Section 4.2 and Wie25 Sections 2.2 and 4 were read). Records were checked on the arXiv abstract page or the DOI page; "full text not consulted" marks a source whose text we did not read.

[AL14] X. An and J. Luk, *Trapped surfaces in vacuum arising dynamically from mild incoming radiation*, [arXiv:1409.6270](https://arxiv.org/abs/1409.6270).*

[BT87] V. B. Braginsky and K. S. Thorne, *Gravitational-wave bursts with memory and experimental prospects*, Nature 327 (1987) 123-125.*

[CCM11] Y. Choquet-Bruhat, P. T. Chruściel and J. M. Martín-García, *The Cauchy problem on a characteristic cone for the Einstein equations in arbitrary dimensions*, Ann. Henri Poincaré 12 (2011) 419-482, [arXiv:1006.4467](https://arxiv.org/abs/1006.4467), [doi:10.1007/s00023-011-0076-5](https://doi.org/10.1007/s00023-011-0076-5) (full text not consulted).

[CFL23] L. Ciambelli, L. Freidel and R. G. Leigh, *Null Raychaudhuri: canonical structure and the dressing time*, [arXiv:2309.03932](https://arxiv.org/abs/2309.03932).*

[Cha26] V. Chandrasekaran, *An algebraic quantization of gravity from the null cone*, [arXiv:2611.34924](https://arxiv.org/abs/2611.34924) (v1, unrefereed; Sections 3 to 5 and Appendices C and D read; equation numbers are those of v1).

[Chr08] D. Christodoulou, *The formation of black holes in general relativity*, [arXiv:0805.3880](https://arxiv.org/abs/0805.3880).*

[Chr91] D. Christodoulou, *Nonlinear nature of gravitation and gravitational-wave experiments*, Phys. Rev. Lett. 67 (1991) 1486-1489.*



[Glo23] J. Glowney, *Arithmoí: A Study of Mathematical Cohesiveness in Physics and Metaphysics*, Book 1 (Volume I, Part 1), First Edition, Athene Noctua, 2023; Sections 1.3.1 (Curvature and Light), 1.3.6 (Hyperbolic Curvature), 1.4.2 and 1.4.4.

[KLR14] S. Klainerman, J. Luk and I. Rodnianski, *A fully anisotropic mechanism for formation of trapped surfaces in vacuum*, Invent. Math. 198 (2014) 1-26, [arXiv:1302.5951](https://arxiv.org/abs/1302.5951), [doi:10.1007/s00222-013-0496-6](https://doi.org/10.1007/s00222-013-0496-6) (full text not consulted).

[KR09] S. Klainerman and I. Rodnianski, *On the formation of trapped surfaces*, [arXiv:0912.5097](https://arxiv.org/abs/0912.5097).*

[Luk12] J. Luk, *On the local existence for the characteristic initial value problem in general relativity*, Int. Math. Res. Not. IMRN 2012, 4625-4678, [arXiv:1107.0898](https://arxiv.org/abs/1107.0898), [doi:10.1093/imrn/rnr201](https://doi.org/10.1093/imrn/rnr201) (full text not consulted).

[Rei07, Rei12] M. P. Reisenberger, *The Poisson bracket on free null initial data for gravity*, [arXiv:0712.2541](https://arxiv.org/abs/0712.2541); *The symplectic 2-form for gravity in terms of free null initial data*, [arXiv:1211.3880](https://arxiv.org/abs/1211.3880).*

[Rend90] A. D. Rendall, *Reduction of the characteristic initial value problem to the Cauchy problem and its applications to the Einstein equations*, Proc. R. Soc. Lond. A 427 (1990) 221-239, [doi:10.1098/rspa.1990.0009](https://doi.org/10.1098/rspa.1990.0009) (full text not consulted).

[Sachs62] R. K. Sachs, *On the characteristic initial value problem in gravitational theory*, J. Math. Phys. 3 (1962) 908-914.*

[Wie24, Wie25] W. Wieland, *Quantum geometry of the null cone*, [arXiv:2401.17491](https://arxiv.org/abs/2401.17491); *Quantum geometry of the light cone: Fock representation and spectrum of radiated power*, [arXiv:2504.10802](https://arxiv.org/abs/2504.10802).*

