---
title: Zeros of the Riemann zeta function
summary: The Hardy Z-function, Gram points and Turing's method for proving the count of zeros, the Riemann–Siegel formula with Arias de Reyna's error bounds, and Platt's method that isolates thousands of zeros at large height with one FFT-based multi-evaluation.
complexity: $O(\sqrt{t})$ operations per evaluation of $Z(t)$ by the Riemann–Siegel formula; Platt's method costs $O(\sqrt{t} + N \log N)$ for $N$ grid points, amortised over thousands of zeros.
order: 55
---

The nontrivial zeros of $\zeta(s)$ with positive imaginary part are
numbered $\rho_n = \beta_n + i\gamma_n$ by increasing $\gamma_n$ (with
multiplicity). FLINT computes them with proof: it finds sign changes of
the real function $Z(t)$, and uses Turing's method to prove that the
number of zeros up to a given height equals the number of sign changes
found, so that every zero in the range is simple, lies on the critical
line, and has a known index. The approach follows the implementation in
mpmath by Arias de Reyna [@Ari2012], with Platt's method [@Pla2017] for
many zeros at large height.

## The Hardy Z-function {#hardy-z}

The Riemann–Siegel theta function
$\theta(t) = \arg\bigl(\pi^{-it/2}\Gamma(\tfrac14 + \tfrac{it}{2})\bigr)$,
defined by continuity with $\theta(0) = 0$, makes

$$
Z(t) = e^{i\theta(t)} \zeta(\tfrac12 + it)
$$

real for real $t$, with $|Z(t)| = |\zeta(\frac12 + it)|$. Every sign
change of $Z$ is a zero of $\zeta$ on the critical line. The number of
zeros with $0 < \gamma \le t$ is

$$
N(t) = \frac{\theta(t)}{\pi} + 1 + S(t), \qquad S(t) = \frac{1}{\pi} \arg \zeta(\tfrac12 + it),
$$

where the argument is defined by continuous variation from $2$ to
$2 + it$ to $\frac12 + it$, and $S(t) = O(\log t)$ with explicit bounds
[@Tru2014]. `acb_dirichlet_hardy_theta` evaluates $\theta(t)$ with the
log-gamma function, and `acb_dirichlet_hardy_z` evaluates $Z(t)$; for
Dirichlet L-functions the same definitions apply with the root number in
$\theta$.

## Gram points and Gram blocks {#gram}

The $n$-th *Gram point* $g_n$ is the solution of $\theta(g_n) = n\pi$
($g_{-1} \approx 9.67$, $g_0 \approx 17.85$), computed by
`acb_dirichlet_gram_point` from an asymptotic formula in terms of the
Lambert W function, refined by Newton's method. Since
$\theta(t) \approx \frac{t}{2}\log\frac{t}{2\pi} - \frac{t}{2} - \frac{\pi}{8}$,
the Gram points are spaced about $2\pi/\log(t/2\pi)$ apart, like the
zeros. A Gram point is *good* if $(-1)^n Z(g_n) > 0$.

*Gram's law* is the observation that $Z$ usually changes sign exactly
once between consecutive Gram points; it locates correctly the first
126 zeros. A *Gram block* is an interval $[g_n, g_{n+k})$ between good
Gram points with no good Gram points inside; *Rosser's rule* says that a
Gram block of length $k$ contains $k$ zeros, and it holds for all Gram
blocks containing the first 13999526 zeros. Both fail infinitely often, but they
are true often enough to guide the search: `acb_dirichlet_isolate_hardy_z_zero`
uses Gram's law up to $n = 126$, Rosser's rule up to $n = 13999526$
(where it is known to hold), and Turing's method beyond.

## Turing's method {#turing}

Counting sign changes gives lower bounds for $N(t)$; Turing's method
gives the matching upper bound from the average of $S(t)$ [@Tur1953].

::: theorem Turing, Lehman, Trudgian {#thm-turing}
For $t_2 > t_1 \ge 168\pi$,

$$
\left| \int_{t_1}^{t_2} S(t)\,dt \right| \le 2.067 + 0.059 \log t_2.
$$
:::

This is the bound of [@Tru2011], which improves those of [@Leh1970] and
[@Bre1979]. If $Z$ has at least the expected number of sign changes in a
run of consecutive Gram blocks after $g_n$ (each block satisfying
Rosser's rule), then $S$ cannot be much larger than $0$ at $g_n$ without
violating the integral bound, which gives $N(g_n) \le n + 1$; a run of
blocks before $g_p$ gives $N(g_p) \ge p + 1$ symmetrically [@Bre1979].
Together with the sign changes counted in $[g_n, g_p)$, this proves
that $N(g_p) - N(g_n) = p - n$ and that all zeros in the range are
simple and on the line.

`acb_dirichlet_turing_method_bound` returns the number of good Gram
blocks needed at height $g_p$: the ceiling of the smaller of
$0.0061 \log^2 g_p + 0.08 \log g_p$ and
$0.0031 \log^2 g_p + 0.11 \log g_p$, derived from the two integral
bounds. The isolation routines scan Gram blocks outward from the target
until enough consecutive blocks satisfying Rosser's rule are found on
each side, subdividing blocks with missing sign changes until all
expected zeros are separated. `acb_dirichlet_zeta_nzeros` computes
$N(t)$ exactly in the same way, and `acb_dirichlet_backlund_s` derives
$S(t)$ from it.

## The Riemann–Siegel formula {#riemann-siegel}

Euler–Maclaurin summation costs $O(t)$ operations at height $t$ (see
[](algorithms/special-functions/zeta#euler-maclaurin)); the
Riemann–Siegel formula needs $O(\sqrt{t})$. With $a = \sqrt{t/2\pi}$,
$N = \lfloor a \rfloor$ and $p = 1 - 2(a - N)$, for $s = \sigma + it$,

$$
\zeta(s) = \mathcal{R}(s) + X(s) \overline{\mathcal{R}}(1-\bar s),
\qquad X(s) = \pi^{s-1/2} \frac{\Gamma((1-s)/2)}{\Gamma(s/2)},
$$

$$
\mathcal{R}(s) = \sum_{k=1}^N \frac{1}{k^s} + (-1)^{N-1} U a^{-\sigma}
\left[ \sum_{k=0}^K \frac{C_k(p)}{a^k} + RS_K \right],
$$

where $U$ is a unimodular factor and the coefficients $C_k(p)$ are
combinations of derivatives of
$F(p) = \cos(\pi(p^2/2 + 3/8))/\cos(\pi p)$ with coefficients
$d_j^{(k)}$ depending on $\sigma$. Arias de Reyna [@Ari2011] gives
recurrences for the $d_j^{(k)}$ (`acb_dirichlet_zeta_rs_d_coeffs`), the
derivatives $F^{(j)}(p)$ by power series division
(`acb_dirichlet_zeta_rs_f_coeffs`), and a rigorous bound for $RS_K$
(Theorem 4.2 of [@Ari2011], `acb_dirichlet_zeta_rs_bound`) valid for
$t$ large enough relative to $K$. The asymptotic series converges only
up to $K = O(t)$ terms, which limits the attainable precision at a given
height; `acb_dirichlet_zeta` uses the formula from
$|t| \ge 24 p^{3/2}$, where it always suffices.

The main cost is the sum of $N \approx \sqrt{t/2\pi}$ powers $k^{-s}$.
It is computed with the sieve of `acb_dirichlet_powsum_sieved`, or, when
the required absolute accuracy (about $p + 10 - \log_2 t$ bits after the
point) fits, with `dfloat_powsum_sieved` in double-word floating-point
arithmetic with rigorous error bounds, which is several times faster.
Derivatives for `acb_dirichlet_zeta_jet_rs` are obtained by a finite
difference.

## Isolating and refining zeros {#refinement}

Given an interval $(a, b)$ in which $Z$ changes sign exactly once,
`_acb_dirichlet_refine_hardy_z_zero` locates the zero. The value of $Z$
near a zero only needs to be known to about $p - \log_2 t$ bits after
the binary point to determine its sign, so $Z$ is evaluated at
$p + 12$ bits while the iterates carry $p + \log_2 t + 8$ bits.

- Below $4\log_2 t + 40$ bits, the Illinois variant of the secant method
  is used, which converges superlinearly while always keeping a
  bracket. When an iterate falls so close to the zero that the sign of
  $Z$ cannot be determined, the zero is bracketed by evaluating at a
  quarter of the tolerance on either side. Within this range, from
  $\log_2 t + 96$ bits on, the zero is first refined in the same way to
  about half the precision, where $Z$ is much cheaper, and then finished
  by one secant step and sign checks at the target precision, four
  evaluations in total.
- From $4 \log_2 t + 40$ bits, Newton's method with
  `acb_dirichlet_hardy_z` of length 2 is used.

The zeros of a range are refined in parallel with `flint_parallel_do`.

## Platt's method {#platt}

At large height, each evaluation of $Z$ by the Riemann–Siegel formula
costs $O(\sqrt{t})$, and isolating $n$ zeros needs several evaluations
per zero. Platt's method [@Pla2011], [@Pla2017], building on the
FFT-based approach of Odlyzko and Schönhage [@OS1988] and its rigorous
version by Booker [@Boo2006], evaluates at once the real function

$$
F(t) = \Lambda(t) e^{\pi t / 4}, \qquad
\Lambda(t) = \pi^{-it/2} \Gamma\left(\tfrac{1/2+it}{2}\right) \zeta\left(\tfrac12 + it\right),
$$

which has the same zeros as $Z(t)$, on a grid of $N = AB$ points
$T - B/2 + k/A$ around a height $T$:

1. The product of $F(T + t)$ with a Gaussian $e^{-t^2/(2h^2)}$ has a
   Fourier transform that is a sum over $j \le J$ of terms
   $j^{-1/2}$ times smooth functions of $x - \log(j\sqrt\pi)/(2\pi)$.
   The terms are grouped by the grid point $m/B$ nearest to
   $\log(j\sqrt\pi)/(2\pi)$ and expanded in Taylor series of length $K$
   around it, which turns the transform into $K$ convolutions of the
   sums $S_{m,k} = \sum_j j^{-1/2} (\log(j\sqrt\pi)/(2\pi) - m/B)^k$ with
   tables of Gaussian-gamma factors.
2. The convolutions are computed with fast Fourier transforms of length
   $N$ (`acb_dft` and [`gr_dft`](gr_dft)), summed in the frequency
   domain, and transformed back, which gives $F$ on the whole grid.
3. Every truncation is bounded rigorously: the Dirichlet series at $J$
   (Lemma B.1 of [@Pla2017]), the Taylor expansion at $K$ (Lemma B.2),
   the aliasing of the discrete transform (Lemmas A.5–A.11), and the
   Gaussian window (Lemma 3.2), with the bounds
   `acb_dirichlet_platt_lemma_*`.
4. Between grid points, $F$ is evaluated by Gaussian-windowed
   Whittaker–Shannon interpolation of the grid values, with rigorous
   bounds for the truncation of the interpolation sum and for the
   window (`acb_dirichlet_platt_ws_interpolation`); this is used to find
   sign changes between grid points and to refine zeros.

The sums $S_{m,k}$ over $j \le J \approx \sqrt{T}$ dominate the cost;
`acb_dirichlet_platt_multieval` computes them in double-word arithmetic
when the precision allows (`acb_dirichlet_platt_use_dfloat`), and
distributes the work over threads. The parameters $A$, $B$, $h$, $J$,
$K$, the interpolation parameters and the odd integers $\sigma$ of the
error bounds are chosen by fitted formulas in terms of $\log T$ for
$10^4 \le n < 10^{37}$. One multi-evaluation isolates thousands of
zeros, and the per-zero cost after it is about 10 ms at $n = 10^{15}$;
`_acb_dirichlet_hardy_z_zeros_use_platt` compares measured costs of both
methods to decide which one `acb_dirichlet_hardy_z_zeros` uses (Platt's
method from $n = 10^{11}$ when enough zeros are requested). The accuracy
of the grid is limited by the parameters, which depend only on the
height, to about 150–250 bits; higher precision is obtained by refining
each zero with the Riemann–Siegel formula.

## The secondary zeta function {#secondary}

$\mathcal{Z}(s) = \sum_{n \ge 1} \gamma_n^{-s}$ converges for
$\operatorname{Re}(s) > 1$ and extends meromorphically to $\mathbb{C}$. Arias de
Reyna [@Ari2020] writes it as $A(s) - P(s) + E(s) - S(s)$, where the
terms are given by convergent series and integrals involving the gamma
function and the prime powers, depending on a free parameter $a > 0$.
`acb_dirichlet_secondary_zeta` implements this decomposition with
rigorous bounds for the tails of $A$, $P$ and $S$ (replacing the
heuristic bounds of [@Ari2020]), without assuming the Riemann
hypothesis.

## Use in FLINT {#use}

- `acb_dirichlet_hardy_theta`, `acb_dirichlet_hardy_z`,
  `acb_dirichlet_gram_point`.
- `acb_dirichlet_zeta_rs`, `acb_dirichlet_zeta_rs_r`,
  `acb_dirichlet_zeta_rs_bound`, `acb_dirichlet_zeta_jet_rs`;
  `dfloat_powsum_sieved` for the main sum.
- `acb_dirichlet_turing_method_bound`, `acb_dirichlet_isolate_hardy_z_zero`,
  `_acb_dirichlet_refine_hardy_z_zero`, `acb_dirichlet_zeta_nzeros`,
  `acb_dirichlet_backlund_s`, `acb_dirichlet_backlund_s_bound`.
- `acb_dirichlet_hardy_z_zeros`, `acb_dirichlet_zeta_zeros` and the
  single-zero versions, which choose between
  `_acb_dirichlet_hardy_z_zeros_rs` and
  `acb_dirichlet_platt_hardy_z_zeros`.
- `acb_dirichlet_platt_multieval`, `acb_dirichlet_platt_ws_interpolation`,
  `acb_dirichlet_platt_local_hardy_z_zeros`.
- `acb_dirichlet_secondary_zeta`.
