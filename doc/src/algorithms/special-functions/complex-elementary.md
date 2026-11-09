---
title: Complex elementary functions
summary: How acb, acf and nfloat_complex compose the real elementary functions into complex ones - numerically stable piecewise formulas, the branch cuts, the per-part rounding of acf with its complex-relative relaxed modes, and the exact, axial and near-axial cases.
complexity: A few evaluations of real elementary functions at the target precision plus a few guard bits; for acf in the strict modes, retries with more guard bits as for arf.
order: 16
---

The complex elementary functions of `acb` (balls), `acf` (rounded
floating point) and `nfloat_complex` (fixed precision) are built from
the real functions of [](algorithms/special-functions/elementary):
each part of the result is a short combination of real exponentials,
logarithms, sines and cosines, chosen so that no part is spoiled by
cancellation in a region where it could be computed accurately. This
page gives the formulas, the branch cuts, and the way `acf` rounds the
parts.

Throughout, $z = a + bi$ with real $a$ and $b$.

## Formulas {#formulas}

**Exponential, sine and cosine.** The parts are products of real
functions of $a$ and of $b$,
$$
e^{a + bi} = e^a(\cos b + i\sin b), \qquad
\sin z = \sin a\cosh b + i\cos a\sinh b, \qquad
\cos z = \cos a\cosh b - i\sin a\sinh b,
$$
and each factor is computed with its own relative accuracy, so that each
part has the relative accuracy of a product. The hyperbolic functions
are $\sinh z = -i\sin(iz)$ and $\cosh z = \cos(iz)$, and the functions
of $\pi z$ use the real functions of $\pi a$ (with the exact reduction of
$a$) and of $\pi b$.

**Logarithm.** $\log z = \log|z| + i\arg z$, with the real part
$\log\sqrt{a^2 + b^2}$ computed by `arb_log_hypot`, which forms
$a^2 + b^2 - 1$ exactly (or to the needed precision) near the unit
circle so that the real part keeps its relative accuracy there, and the
imaginary part by $\operatorname{atan2}(b, a)$.

**Tangent and cotangent.** In `acb`, for $|b| < 1$ the tangent is the
quotient $\sin z/\cos z$, and for $|b| \ge 1$
$$
\tan z = -i\,\frac{e^{2iz} - 1}{e^{2iz} + 1} \quad (b > 0), \qquad
\tan z = i\,\frac{e^{-2iz} - 1}{e^{-2iz} + 1} \quad (b < 0),
$$
where the exponential is small, so that the result $\pm i(1 - 2e^{\pm 2iz}
+ \cdots)$ comes out without cancellation; the cotangent likewise. The
rounded functions of `acf` use the forms
$$
\tan z = \frac{\sin a\cos a + i\sinh b\cosh b}{\cos^2 a + \sinh^2 b},
\qquad
\cot z = \frac{\sin a\cos a - i\sinh b\cosh b}{\sin^2 a + \sinh^2 b},
$$
and $\sec z = (\cos a\cosh b + i \sin a\sinh b)/(\cos^2 a + \sinh^2 b)$,
$\csc z = (\sin a\cosh b - i\cos a\sinh b)/(\sin^2 a + \sinh^2 b)$, whose
denominators $|\cos z|^2$ and $|\sin z|^2$ are sums of squares: each
part is a quotient of accurately known numbers, also when it is
exponentially small against the other (the real part of $\tan z$ is about
$4\sin a\cos a\,e^{-2|b|}$ for large $|b|$). The hyperbolic tangent and its
relatives are $\tanh z = -i\tan(iz)$, $\coth z = -i\cot(-iz)$,
$\operatorname{sech} z = \sec(iz)$ and $\operatorname{csch} z =
-i\csc(-iz)$.

**Square root.** For $a \ge 0$, $\sqrt{z} = \sqrt{(|z| + a)/2} +
i\,b/\sqrt{2(|z| + a)}$, and for $a < 0$, $\sqrt{z} = |b|/\sqrt{2(|z| -
a)} + i\operatorname{sgn}(b)\sqrt{(|z| - a)/2}$: the sum $|z| \pm a$ is
always one without cancellation.

**Inverse functions.** In terms of the logarithm and the square root,
$$
\operatorname{atan} z = \frac{\log(1 + iz) - \log(1 - iz)}{2i}, \qquad
\operatorname{asin} z = -i\log\bigl(iz + \sqrt{1 - z^2}\bigr), \qquad
\operatorname{acos} z = \frac\pi2 - \operatorname{asin} z,
$$
$$
\operatorname{acosh} z = \log\bigl(z + \sqrt{z + 1}\sqrt{z - 1}\bigr),
\qquad \operatorname{asinh} z = -i\operatorname{asin}(iz), \qquad
\operatorname{atanh} z = -i\operatorname{atan}(iz),
$$
with `log1p` for the arctangent so that it stays accurate near 0, and
the real functions on the real axis inside their domains. These
formulas also define the branch cuts below.

On the axes all functions reduce to real ones (for example $\log(-x) =
\log x + \pi i$ and $\log(ix) = \log|x| \pm \pi i/2$ for $x > 0$, and
$\tan(ib) = i\tanh b$), and the parts that are exactly zero are set to
zero.

## Branch cuts {#branch-cuts}

The multivalued functions return their principal values, with the branch
cuts of the DLMF [@NIST2012]: $(-\infty, 0]$ for the logarithm and the
square root, $(-\infty, -1) \cup (1, +\infty)$ for $\operatorname{asin}$,
$\operatorname{acos}$ and $\operatorname{atanh}$, $(-i\infty, -i) \cup
(i, +i\infty)$ for $\operatorname{atan}$ and $\operatorname{asinh}$, and
$(-\infty, 1)$ for $\operatorname{acosh}$; the values on a cut are those
of [](acb/branch-cuts#conventions). In ball arithmetic a rectangle that
touches or crosses a cut must enclose the values on both sides, so its
result is wide; the analytic variants of `acb` return an indeterminate
result instead ([](acb/branch-cuts)). `acf` has no signed zeros: a
number with a zero imaginary part lies on the cut and gets the value
given there.

## Register paths {#fast}

On 64-bit machines, at precisions of one and two limbs, `acb_exp`,
`acb_exp_pi_i`, `acb_log`, `acb_sin_cos` and `acb_sin_cos_pi` of
narrow rectangles, the functions of `nfloat_complex` and those of `acf`
with mantissas of at most two limbs compose the register paths of the
real kernels ([](algorithms/special-functions/elementary#fast)) into
both parts at once, with an error count per part. For `acb` the input
radii are then propagated through bounds on the derivative, as for the
real functions; for `acf` the parts are certified like the frames of
`arf`.

## Rounding in acf {#acf}

`acf` rounds each part of the result to `prec` bits. In the strict modes
each part is correctly rounded, whatever its size against the other:
the compositions are evaluated at an increasing working precision until
the error bound of each part decides its rounding, as in
[](algorithms/floating-point/correct-rounding). Since every part is a
product or quotient of accurately known real numbers (or a real function
of an exactly known number), its relative error bound is a small
multiple of the working precision, independently of the other part.

The relaxed modes guarantee the relative accuracy of the complex number
instead: in `ARF_RND_FAST` each part is within 1 ulp and in
`ARF_RND_ACCURATE` within 0.51 ulp of the larger part (at `prec` bits).
A part much smaller than the other may then have a worse relative
accuracy or be zero, but such a part does not slow the evaluation down.
In the strict modes, a part exponentially smaller than the other must be
resolved; for the tangent and cotangent near $|\operatorname{Im}| = 1$,
the imaginary part is $\pm(1 + \delta)$ with
$$
\delta = -\frac{\cos 2a + e^{-2|b|}}{2(\cos^2 a + \sinh^2 b)}
\quad\text{resp.}\quad
\frac{\cos 2a - e^{-2|b|}}{2(\sin^2 a + \sinh^2 b)},
$$
which is computed directly, so that the directed modes are rounded from
the sign of $\delta$ when $|\delta|$ is far below an ulp: $\cot(2^k i -
1)$ takes about a microsecond at 53 bits in every mode for $k$ up to
about 1000 (beyond, the cost grows with the bits of $e^{2|b|}$ that the
real part needs: a few milliseconds at $k = 65536$). The tangent family
gives nan for $|\operatorname{Im} z|$ (the real part for the hyperbolic
functions) beyond about $2^{131072}$, like the exponential functions it
is made of.

**Exact results.** No window of bits decides the rounding of a dyadic
value, so the exact results must be recognised: the zeros of the sine
and cosine of $\pi x$ at the half-integers (the imaginary part of
$(-2)^3$, say), integer powers (and small dyadic powers of exact roots)
with at most $2^{24}$ bits, the powers of $2^j(\pm 1 \pm i)$ for any
integer exponent, dyadic square roots, and the parts that are exactly
zero or a multiple of $\pi/2$ on the axes.

**Values next to a rounding boundary.** A part can also lie within any
distance of a dyadic number without being one, when the parts of the
argument differ in size by many bits: the real part of $1/(1 +
2^{-k}i)$ is $1/(1 + 2^{-2k})$, that of $e^{\varepsilon(1 + i)}$ is $1 +
\varepsilon + O(\varepsilon^3)$. Retrying until the gap shows would take a
working precision proportional to it. Instead, a part still undecided
at a working precision well above `prec` is compared exactly with the
rounding boundary $\beta$ that its approximation is next to (the
nearest number of `prec` bits, or of `prec` + 1 bits for `ARF_RND_NEAR`):
the part rounds like $\beta \pm$ a sliver for the sign of $V - \beta$, or
like $\beta$ for $V = \beta$.

- For `acf_inv`, `acf_div`, `acf_sqrt` and `acf_rsqrt` the sign of $V -
  \beta$ is that of a polynomial in the inputs and $\beta$: $a - \beta(a^2
  + b^2)$ for the real part of $1/(a + bi)$; for the square root $q =
  \sqrt{(r + |a|)/2}$, $q > \beta \Leftrightarrow r > 2\beta^2 - |a|$, which
  squares to $b^2 - 4\beta^4 + 4\beta^2|a| > 0$ (and likewise for $|b|/(2q)$
  and the quotients by $r$ of the reciprocal root). These are sums of
  exact terms, whose sign is found by adding the terms in decreasing
  magnitude exactly until the rest cannot reach half the partial sum, so
  that no sum spans the gap to a negligible term. Such a part is decided
  at `prec` + 64 bits whatever the gap.
- For the elementary functions, a transcendental part can come that
  close to a dyadic number only near an exact point $z_0$ of the function
  (0; $\pm 1, \pm i$ for the logarithm; the multiples of $1/4$ for the
  functions of $\pi z$; $\pm 1$ for the square-root branch points of the
  inverse sine and cosine; the simple poles of cot, csc and their
  relatives; infinity for the inverse functions), where
  $$
  f(z_0 + h) = C + u\,s \sum_{j = j_0}^{N} q_j (\pi^k H)^j + O(|H|^{N+1})
  $$
  with $H$ a rotation of $h$ (or $1/z$ at infinity), rational $q_j$ (the
  series of `fmpq_poly`), $j_0 = -1$ at a pole, $C \in \mathbb{Q}(i)\pi^c$,
  a unit $u$ and $s = \sqrt{w h}$ at a branch point when it is dyadic.
  Multiplied by a positive denominator ($\pi^k|H|^2$ at a pole,
  $|z|^{2N}$ at infinity), $\operatorname{Re}$ or $\operatorname{Im}$ of
  the expansion minus $\beta$ is a sum of exact monomials in the parts of
  $H$ and $s$ with integer coefficients and powers of $\pi$: summed
  with the exact sign within each power of $\pi$, at any spread of the
  exponents (which catches exact cancellations such as
  $e^{2^{-2k-1} + 2^{-k}i}$, whose real part is $1 - 2^{-4k}/12 + \ldots$),
  correctly rounded by `arf_sum` where several powers of $\pi$ combine,
  and with $\pi$ as a ball across them, the sign decided when the
  sum exceeds the bound on the omitted terms. The imaginary part of the
  logarithm is $\operatorname{atan2}(b, a)$, rounded by `arf` directly.
  The power $x^y$ near $x = 2^K i^u(1 + h)$ is $x_0^y\sum_j \binom{y}{j}
  h^j$ with the binomial sum as exact terms and $x_0^y$ exact when it is
  dyadic, and for a tiny $y$ it is $1 + \operatorname{expm1}(y \log x)$ in
  ball arithmetic.

A value that neither decides (none is known) would keep the retries from
converging, so they stop at a working precision of $4\,\mathit{prec} +
65536$ bits (plus twice the spread of the exponents of the argument's
parts, up to $2^{20}$ bits, for the functions computed from balls
below), and the approximation is then rounded as it stands.

**Inverse functions.** The inverse trigonometric and hyperbolic
functions, `expm1`, `log1p` and the $\operatorname{sinc}$ functions round
the balls of the corresponding `acb` functions at increasing working
precisions: in the strict modes until both endpoints of the ball of each
part round to the same number, in the relaxed modes until the radii lie
below $2^{-\mathit{prec} - 4}$ resp. $2^{-\mathit{prec} - 8}$ times the
larger part. The precision grows at least geometrically, so that a part
that vanishes to many more bits than the accuracy shows is resolved in
a logarithmic number of steps. For $e^z - 1$ with $\operatorname{Re} z
\le -(\mathit{prec} + 3)$ the real part is $-(1 + \delta)$ with a tiny
$\delta = -e^{\operatorname{Re} z}\cos(\operatorname{Im} z)$, rounded from
the sign of $\delta$ like the asymptotic cases of `arf`. Off the axes the parts are transcendental at dyadic
arguments, so that the strict loop terminates. An argument very close to
an axis ($z = z_0 + \eta$ with $|\eta| < 2^{-\mathit{prec} - 64}|z_0|$,
$z_0$ on the axis and away from the branch cuts) is evaluated as
$f(z_0) + f'(z_0)\eta$, with the remainder bounded by the Cauchy
estimates on a disk around $z_0$
([](algorithms/ball-arithmetic/error-propagation#cor-cauchy)), so that
the small part comes out with its own relative accuracy at the usual
working precision: $\operatorname{asin}(1/2 + 2^{-100000} i)$ takes
about two microseconds at 53 bits.

**Division and square roots.** $1/z$ and $z/w$ form the numerators
$a c + b d$, $b c - a d$ and the denominator $c^2 + d^2$ exactly, so
that each part is one correctly rounded quotient of exact numbers (when
the terms of a sum are so far apart that they cannot cancel, they are
rounded to a working precision instead and the quotient is certified).
The square root uses the formulas above at an increasing working
precision, with the dyadic roots computed exactly.

## Fixed precision: nfloat_complex {#nfloat}

The functions of `nfloat_complex` evaluate the same kind of formulas
once, in the precision of the context, with no guard bits except for the
exact evaluation of $|z|^2 - 1$ near the unit circle: each part
typically has an error of a few ulps relative to the magnitude of the
result (a small part of a large result can have a large relative
error), and directed rounding is not respected.

## Use in FLINT {#use}

::: functions acb_exp acb_log acb_sin_cos acb_tan acb_atan acb_asin acb_sqrt acf_exp acf_log acf_tan acf_asin acf_pow acf_sqrt :::
