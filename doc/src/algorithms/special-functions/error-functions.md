---
title: Error functions, incomplete gamma functions and exponential integrals
summary: How the error functions, Fresnel integrals, incomplete gamma and beta functions and the exponential, logarithmic, sine and cosine integrals are reduced to 1F1, U and 2F1, how FLINT chooses between their convergent and asymptotic expansions, the real-variable algorithms with rectangular splitting, and the bit-burst algorithm for the incomplete gamma function at very high precision.
order: 44
---

The error function, the incomplete gamma function and the exponential
integrals are confluent hypergeometric functions with one parameter
fixed or tied to the argument:
$$
\gamma(s, z) = \frac{z^s}{s}\, {}_1F_1(s; s+1; -z), \qquad
\Gamma(s, z) = e^{-z} U(1-s, 1-s, z) = z^s e^{-z} U(1, 1+s, z),
$$
$$
\operatorname{erf}(z) = \frac{\gamma(\tfrac12, z^2)}{\sqrt\pi}\,\frac{z}{\sqrt{z^2}},
\qquad E_s(z) = z^{s-1}\Gamma(1-s,z), \qquad
\operatorname{Ei}(z) = -E_1(-z) + \text{(branch terms)}.
$$
Their evaluation therefore follows
[](algorithms/special-functions/hypergeometric): a convergent series
near the origin, the asymptotic expansion of $U$ with Olver's error
bound for large $|z|$, and a working precision chosen from the
cancellation of the series.

## Error functions {#erf}

FLINT uses three representations of the error function:
$$
\operatorname{erf}(z) = \frac{2z}{\sqrt{\pi}}\,
{}_1F_1(\tfrac{1}{2}, \tfrac{3}{2}, -z^2)
= \frac{2z e^{-z^2}}{\sqrt{\pi}}\, {}_1F_1(1, \tfrac{3}{2}, z^2)
= \frac{z}{\sqrt{z^2}} - \frac{e^{-z^2}}{z \sqrt{\pi}}
U^{*}(\tfrac{1}{2}, \tfrac{1}{2}, z^2).
$$
For $z = x + yi$ the terms of the first series have size up to about
$e^{x^2}$ when $|x| > |y|$ and the function has size about $1$, so the
series loses about $x^2/\log 2$ bits; the second series has the
opposite behaviour (by Kummer's transformation), losing about
$y^2/\log 2$. `acb_hypgeom_erf` uses the second series if $|x| \ge |y|$
and the first otherwise, at a working precision raised by $2
\min(x^2, y^2)/\log 2$ bits, and the asymptotic expansion when its
error term $\approx e^{-|z|^2}/|z|$ is below $2^{-p}$, computing the
$U^*$ term at the lower precision that its small prefactor allows. For
inexact $z$, the series is evaluated at the midpoint and
`acb_hypgeom_erf_propagated_error` bounds the propagated error with the
derivative $\frac{2}{\sqrt\pi} e^{-z^2}$. The complementary function
$\operatorname{erfc}(z) = 1 - \operatorname{erf}(z)$ is computed from
the $U^*$ representation without the subtraction when
$\operatorname{Re}(z)$ is large and positive, and
$\operatorname{erfi}(z) = -i\operatorname{erf}(iz)$. The normalized
Fresnel integrals satisfy [@NIST2012, 7.5]
$$
C(z) \pm i S(z) = \frac{1 \pm i}{2} \operatorname{erf}\left(\frac{\sqrt\pi}{2}(1 \mp i) z\right),
$$
and `acb_hypgeom_fresnel` computes them from $\operatorname{erf}((1 \pm
i) w)$ for a real multiple $w$ of $z$, so that the square $\pm 2iw^2$
inside the error function is computed without rounding of the real
part.

### Real arguments {#erf-real}

For real $x$, `arb_hypgeom_erf` uses
$$
\operatorname{erf}(x) = \frac{2x}{\sqrt\pi} e^{-x^2} \sum_{k=0}^\infty
\frac{(x^2)^k}{(3/2)_k},
$$
whose terms are positive, so that no precision is lost. Since
$(x^2)^k/(3/2)_k \le (x^2)^k/k!$, the tail after $N$ terms is at most
the tail of the series of $e^{x^2}$, and $N$ is chosen with the Lambert
$W$ function. The sum is computed by rectangular splitting
(`_arb_hypgeom_gamma_lower_sum_rs_1`): with $m \approx \sqrt N$ powers
of $x^2$, the coefficients $1/(3/2)_k$ are accumulated as products of
the word-size factors $2/(2k+1)$, and only about $2\sqrt N$
multiplications of full-precision numbers are needed. For $|x| > 2$
the asymptotic expansion
$$
\operatorname{erfc}(x) = \frac{e^{-x^2}}{x\sqrt\pi} \left(\sum_{k=0}^{N-1}
\frac{(1/2)_k}{(-x^2)^k} + R_N\right), \qquad |R_N| \le
\frac{(1/2)_N}{x^{2N}},
$$
whose remainder is bounded by the first omitted term for real $x$, is
used when it reaches the target accuracy; its sum is also computed by
rectangular splitting (`_arb_hypgeom_gamma_upper_sum_rs_1`).
`arb_hypgeom_erfc` uses the asymptotic expansion whenever possible and
otherwise the convergent series at a working precision increased by
$(x^2 + \log x)/\log 2$ bits.

## Incomplete gamma functions {#incomplete-gamma}

The upper incomplete gamma function is computed from
$$
\Gamma(s,z) = e^{-z} U(1-s,1-s,z) = \Gamma(s) - \frac{z^s}{s}\,
{}_1F_1(s, s+1, -z) = \Gamma(s) - \frac{z^s e^{-z}}{s}\,{}_1F_1(1, s+1, z).
$$
The first form is the asymptotic expansion
$\Gamma(s,z) \sim z^{s-1}e^{-z}\sum_k (1-s)_k(-z)^{-k}$, whose error
bound from [](algorithms/special-functions/hypergeometric#asymptotic)
depends on $r = |b - 2a| = |s - 1|$. The two convergent forms differ by
Kummer's transformation; the second is used for $\operatorname{Re} z
\le 0$ and the third for $\operatorname{Re} z > 0$, so that the series
has terms of equal sign on the real axis. For $s = -n$ a nonpositive
integer, $\Gamma(s)$ has a pole that cancels with the series, and the
limit
$$
\Gamma(-n,z) = \frac{(-1)^n}{n!} (\psi(n+1) - \log(z))
+ \frac{(-1)^n}{(n+1)!} z \, {}_2F_2(1,1;2,2+n;-z)
- z^{-n} \sum_{k=0}^{n-1} \frac{(-z)^k}{(k-n) k!}
$$
is used instead. `acb_hypgeom_gamma_upper` chooses the asymptotic
expansion for nonnegative real $s$ and $z$ when $(z - 2.39)^8 > (1.03
s)^8 + (0.33\,p)^8$ (an empirical fit of the region where it reaches
$p$ bits), and for other arguments when `acb_hypgeom_u_use_asymp`
holds and $s$ is a positive integer or $z$ lies in one of Olver's
regions. The lower incomplete gamma function is
$\gamma(s,z) = \frac{z^s}{s}{}_1F_1(s; s+1; -z)$, computed with
`acb_hypgeom_m` (regularized when requested), which applies the same
choices. The regularized forms $Q(s,z) = \Gamma(s,z)/\Gamma(s)$,
$P(s,z) = \gamma(s,z)/\Gamma(s)$ and $\gamma^*(s,z) = z^{-s}P(s,z)$ are
computed directly with regularized hypergeometric functions, so that
they are defined at the poles of $\Gamma(s)$.

### Rational parameters and the bit-burst algorithm {#bit-burst}

For a rational parameter $a = p/q$ and real $z > 0$, the internal
functions `_arb_hypgeom_gamma_upper_fmpq_inf_bsplit` and
`_arb_hypgeom_gamma_lower_fmpq_0_bsplit` sum the asymptotic expansion
at infinity and the Taylor series
$$
\gamma(a, z) = \frac{z^a e^{-z}}{a} \sum_{n=0}^\infty
\frac{z^n}{(a+1)_n}
$$
at zero by binary splitting, with truncation errors bounded by the
corresponding `_choose_N` functions for a prescribed absolute
tolerance. This is efficient only when $z$ has few bits. For a
full-precision $z$, `arb_hypgeom_erf_bb` applies the bit-burst
algorithm ([](algorithms/special-functions/hypergeometric#bit-burst))
to $\Gamma(1/2, z^2)$:

1. $z_0$ is $z^2$ truncated to 64 fractional bits (or more if $z^2$ is
   tiny), and $\Gamma(a, z_0)$ is computed with the asymptotic series
   if it reaches the tolerance and as $\Gamma(a) - \gamma(a, z_0)$
   otherwise.
2. For $j = 1, 2, \ldots$, $z_j$ is $z^2$ truncated to $2^{j} \cdot 64$
   bits, and `_arb_gamma_upper_fmpq_step_bsplit` computes $\Gamma(a,
   z_j)$ from $\Gamma(a, z_{j-1})$ and $e^{-z_{j-1}}$ with the Taylor
   expansion at $z_{j-1}$; the coefficients satisfy a recurrence from
   the differential equation $z y'' + (1 + z - a) y' = 0$ of $y =
   \Gamma(a, z)$, so the step is a product of $2 \times 2$ polynomial
   matrices computed by binary splitting, with a Cauchy-type bound for
   the truncation error.
3. A last step goes from the last truncation to $z^2$ itself.

The working precisions of the two parts are adjusted to the absolute
tolerance needed for $\operatorname{erf}$ or $\operatorname{erfc}$
(with cancellation when $\operatorname{erfc}$ is computed from the
lower series). `arb_hypgeom_erf` and `arb_hypgeom_erfc` use this
method above 30000 bits for moderately large $|x|$; it costs
$O(M(p)\log^2 p)$, compared with $O(\sqrt{p})$ full multiplications
plus $O(p)$ cheap ones for rectangular splitting. The same functions
are used by [`acb_dirichlet`](acb_dirichlet) for the incomplete gamma
functions in the approximate functional equation of Dirichlet
$L$-functions at rational $s$.

## Incomplete beta function {#beta}

$$
B(a,b;z) = \frac{z^a}{a}\, {}_2F_1(a, 1-b, a+1, z), \qquad
I(a,b;z) = \frac{\Gamma(a+b)}{\Gamma(b)}\, z^a\,
{}_2\mathbf{F}_1(a, 1-b, a+1, z)
$$
are computed with `acb_hypgeom_2f1`, and inherit its transformations
and analytic continuation ([](algorithms/special-functions/gauss-2f1)).

## Exponential, logarithmic and trigonometric integrals {#expint}

The generalized exponential integral is $E_s(z) = z^{s-1}\Gamma(1-s,
z)$ (`acb_hypgeom_expint`). The other integrals have the
representations
$$
\operatorname{Ei}(z) = z\, {}_2F_2(1, 1; 2, 2; z) + \gamma +
\tfrac{1}{2} \left(\log(z) - \log(1/z)\right)
= -e^z U(1,1,-z) - \log(-z) + \tfrac{1}{2} \left(\log(z) - \log(1/z)\right),
$$
$$
\operatorname{Si}(z) = z\, {}_1F_2(\tfrac{1}{2}; \tfrac{3}{2},
\tfrac{3}{2}; -\tfrac{z^2}{4}) = \frac{i}{2} \left[ e^{iz} U(1,1,-iz) -
e^{-iz} U(1,1,iz) + \log(-iz) - \log(iz) \right],
$$
$$
\operatorname{Ci}(z) = \log(z) + \gamma - \tfrac{z^2}{4}\,
{}_2F_3(1, 1; 2, 2, \tfrac{3}{2}; -\tfrac{z^2}{4}) = \log(z) -
\tfrac{1}{2} \left[ e^{iz} U(1,1,-iz) + e^{-iz} U(1,1,iz) + \log(-iz) +
\log(iz) \right],
$$
and similar ones for $\operatorname{Shi}(z) = -i\operatorname{Si}(iz)$
and $\operatorname{Chi}$, while $\operatorname{li}(z) =
\operatorname{Ei}(\log z)$. The logarithmic terms implement the branch
cuts of *Mathematica*: $\operatorname{Ei}$, $\operatorname{Ci}$ and
$\operatorname{Chi}$ have a cut on $(-\infty, 0]$, and the combination
$\frac12(\log z - \log(1/z))$ equals $\log z$ except on the cut, where
it averages the two sides. For each function the automatic version
compares the accuracy of the asymptotic expansion ($|z|/\log 2$ bits)
with that of the series after cancellation, as for $U$ and $M$.

For real arguments, `arb_hypgeom_si` and `arb_hypgeom_ci` avoid complex
arithmetic: the convergent series are summed with
`arb_hypgeom_sum_fmpq_arb` at a working precision raised by $|z|/\log
2$ bits, and the asymptotic expansions
$$
\operatorname{Si}(z) = \frac\pi2 - \frac{\cos z}{z} f(z) - \frac{\sin z}{z}
g(z), \qquad \operatorname{Ci}(z) = \frac{\sin z}{z} f(z) - \frac{\cos
z}{z} g(z),
$$
where $f(z) + i g(z) \sim \sum_k k!\,(i/z)^{k}$, with `arb_hypgeom_sum_fmpq_imag_arb`
in powers of $1/z$; the error of the asymptotic expansion is bounded by
the first omitted term $N!/z^N$. The number of terms of the asymptotic
expansion is the first $N$ for which this bound is below $2^{-p}$, if
the bound reaches it before it starts to grow.

## Use in FLINT {#use}

- Error functions: `acb_hypgeom_erf`, `acb_hypgeom_erfc`,
  `acb_hypgeom_erfi` and the methods `acb_hypgeom_erf_1f1a`,
  `acb_hypgeom_erf_1f1b`, `acb_hypgeom_erf_asymp`; real versions
  `arb_hypgeom_erf`, `arb_hypgeom_erfc` with `arb_hypgeom_erf_bb`;
  inverses `arb_hypgeom_erfinv`, `arb_hypgeom_erfcinv` (Newton
  iteration, see [](arb_hypgeom/erf)); Fresnel integrals
  `acb_hypgeom_fresnel`.
- Incomplete gamma and beta functions: `acb_hypgeom_gamma_upper` and
  its methods, `acb_hypgeom_gamma_lower`, `acb_hypgeom_beta_lower`,
  the real wrappers, `arb_hypgeom_gamma_upper_integration`, and the
  rational-parameter functions `_arb_hypgeom_gamma_upper_fmpq_inf_bsplit`,
  `_arb_hypgeom_gamma_lower_fmpq_0_bsplit`,
  `_arb_gamma_upper_fmpq_step_bsplit`.
- Exponential integrals: `acb_hypgeom_expint`, `acb_hypgeom_ei`,
  `acb_hypgeom_si`, `acb_hypgeom_ci`, `acb_hypgeom_shi`,
  `acb_hypgeom_chi`, `acb_hypgeom_li` and the real versions, with
  `arb_hypgeom_si` and `arb_hypgeom_ci` implemented directly.
- All of them have power series versions (`_series`), which compute
  the constant term with the scalar function and the other
  coefficients by integrating the derivative, which is elementary:
  $\frac{d}{dz}\operatorname{erf}(z) = \frac{2}{\sqrt\pi}e^{-z^2}$,
  $\frac{d}{dz}\Gamma(s,z) = -z^{s-1}e^{-z}$, $\frac{d}{dz}
  \operatorname{Si}(z) = \sin(z)/z$, and so on.
