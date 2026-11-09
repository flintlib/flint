---
title: acf.h
header: acf.h
description: Complex floating-point numbers as pairs of arf_t, for approximate computation
category: Real and complex numbers
order: 5
---

An `acf_t` is a complex floating-point number $a + bi$ whose real and
imaginary parts are `arf_t` numbers, without error bounds. It is the
complex counterpart of `arf_t` and is used where complex arithmetic is
needed but error bounds are not: in approximate (non-rigorous)
algorithms such as the floating-point stage of eigenvalue or root
computations that are certified afterwards, and as the element type of the
generic ring of complex floating-point numbers
(`gr_ctx_init_complex_float_acf`).

Addition, subtraction and multiplication are
correctly rounded in each part, like the `arf` operations; division,
square root and dot products are only approximate. All functions allow
aliasing of inputs and outputs.

## Types {#types}

::: type
acf_struct
acf_t
acf_ptr
acf_srcptr

An `acf_struct` has two `arf_struct` fields, the real and imaginary
parts. `acf_t` is an array of one `acf_struct`, passed by reference;
`acf_ptr` and `acf_srcptr` point to arrays.
:::

::: macro
acf_realref(x)
acf_imagref(x)

Pointers to the real and imaginary parts of `x` as `arf_t`.
:::

::: function
arf_ptr acf_real_ptr(acf_t z)
arf_ptr acf_imag_ptr(acf_t z)

Function versions of `acf_realref` and `acf_imagref`, for language
bindings.
:::

## Memory management {#memory}

::: function
void acf_init(acf_t x)
void acf_clear(acf_t x)

Initialise `x` to zero, and free the memory used by `x`.
:::

::: function
acf_ptr _acf_vec_init(slong n)
void _acf_vec_clear(acf_ptr v, slong n)

Allocate a vector of `n` numbers initialised to zero, and free such a
vector.
:::

::: function
void acf_swap(acf_t z, acf_t x)

Swaps `z` and `x` without copying.
:::

::: function
slong acf_allocated_bytes(const acf_t x)

Returns the number of bytes allocated on the heap by `x`, excluding
`sizeof(acf_struct)`.
:::

## Basic operations {#basic}

::: function
void acf_set(acf_t z, const acf_t x)
void acf_neg(acf_t z, const acf_t x)

Set `z` to `x`, respectively $-x$, exactly.
:::

::: function
int acf_set_round(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_neg_round(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)

Set `res` to `x`, respectively $-x$, with both parts rounded to `prec`
bits in the direction `rnd`. Bit 0 of the return value is the inexact
flag of the real part and bit 1 that of the imaginary part.
:::

::: function
int acf_equal(const acf_t x, const acf_t y)

Returns whether `x` and `y` are equal (both parts equal in the sense of
`arf_equal`).
:::

::: function
slong acf_bits(const acf_t x)

Returns the larger of `arf_bits` of the two parts.
:::

::: function
void acf_get_mag(mag_t res, const acf_t x)

Sets `res` to an upper bound for $|x|$.
:::

::: function
void acf_printd(const acf_t x, slong n)

Prints `x` as `a + b*I` with `n` significant decimal digits for each part
(see `arf_printd`).
:::

::: function
void acf_randtest(acf_t x, flint_rand_t state, slong bits, slong mag_bits)

Sets both parts of `x` to random numbers as `arf_randtest` does.
:::

## Arithmetic {#arithmetic}

::: function
int acf_add(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
int acf_sub(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
int acf_mul(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)

Set `res` to $x + y$, $x - y$, respectively $xy$, with each part correctly
rounded to `prec` bits in the direction `rnd`. Bit 0 of the return value
is the inexact flag of the real part and bit 1 that of the imaginary part.

::: details
`acf_mul` uses `arf_complex_mul`, or `arf_complex_sqr` when `x` and `y`
are the same object; each part of the result is the exact $ac - bd$ or
$ad + bc$ rounded once.
:::
:::

## Approximate arithmetic {#approx}

The following functions are not correctly rounded: `rnd` is the direction
of the final rounding, but the intermediate roundings are not specified
and no error bound is given.

::: function
void acf_approx_inv(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
void acf_approx_div(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)

Set `res` to an approximation of $1/x$, respectively $x/y$. A real or
purely imaginary divisor is handled by real divisions of the parts;
otherwise, above two limbs of precision, the complex kernels of
[](arf#mpn) compute the parts within about an ulp (relative to the
magnitude of the result for the quotient, whose numerator parts may
cancel) at $\lceil \mathit{prec} / B \rceil$ limbs, rounded to `prec`
bits in the direction `rnd`, for operands with moderate exponents;
otherwise $1/(a + bi) = (a - bi)/(a^2 + b^2)$ with $a^2 + b^2$ computed
by `arf_sosq`.
:::

::: function
void acf_approx_sqrt(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)

Sets `res` to an approximation of the principal square root of `x`,
through the complex kernel of [](arf#mpn) for operands with
nonzero parts and moderate exponents.
:::

::: function
void acf_approx_rsqrt(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)

Sets `res` to an approximation of $1 / \sqrt{x}$ (the principal root),
through the complex kernel of [](arf#mpn) for operands with nonzero
parts and moderate exponents, otherwise as the reciprocal of the square
root; `rsqrt(0)` is $+\infty$.
:::

::: function
void acf_approx_dot(acf_t res, const acf_t initial, int subtract, acf_srcptr x, slong xstep, acf_srcptr y, slong ystep, slong len, slong prec, arf_rnd_t rnd)

Sets `res` to an approximation of $s + (-1)^{\mathit{subtract}}
\sum_{i=0}^{\mathit{len}-1} x_i y_i$, where $s$ is `initial` (0 if
`NULL`), with the conventions of `arb_dot`, using a fixed-point
accumulation like that of `acb_approx_dot`.
:::

## Rounded arithmetic and elementary functions {#elem}

The following functions round each part of the result to `prec` bits.
In the strict modes (`ARF_RND_DOWN`, `ARF_RND_UP`, `ARF_RND_FLOOR`,
`ARF_RND_CEIL`, `ARF_RND_NEAR`) each part is correctly rounded, as the
functions of [](arf#elementary) round their results, whatever the size
of the part against the other one. The relaxed modes guarantee the
relative accuracy of the complex number instead: in `ARF_RND_FAST` each
part is within 1 ulp and in `ARF_RND_ACCURATE` within 0.51 ulp of the
larger part (at `prec` bits) of the exact value, so that a part much
smaller than the other may have a worse relative accuracy or be zero.
They are not slowed down by a part that is exponentially small against
the other one, such as the real part of $\cot(2^k i - 1)$, which the
strict modes must resolve. The return value carries the inexact flags
of the parts in the strict modes (bit 0 the real part, bit 1 the
imaginary part; bits 2 and 3 those of a second result); in the relaxed
modes it is not a reliable inexact flag.

The parts are computed from the real functions of [](arf#elementary)
by formulas that keep each part accurate, including a part that is
exponentially smaller than the other, at an increasing working
precision until the rounding of each part is decided (in the strict
modes). Exact results are recognised (for example the zeros of
$\sin(\pi x)$, exact integer powers and dyadic square roots), and a
part still undecided at a working precision well above `prec` (one that
lies extremely close to a rounding boundary, as when the parts of the
argument differ in size by very many bits) is decided exactly from its
side of that boundary: by a polynomial identity for `acf_inv`,
`acf_div`, `acf_sqrt` and `acf_rsqrt`, and from the expansion of the
function at its nearby exact point (0, $\pm 1$, $\pm i$, the multiples of
1/4 for the functions of $\pi z$, the branch points and poles, infinity)
for the others, without a precision growing with that distance. A value
that neither decides (none is known) would keep the retries from
converging, so beyond a working precision of $4 \mathit{prec} + 65536$
bits (plus up to $2^{21}$ bits for arguments whose parts have very
different exponents) the approximation is rounded as it stands. The tangent family gives nan beyond about
$|\operatorname{Im} z| = 2^{131072}$. The formulas, the exact cases and the rounding are described
in [](algorithms/special-functions/complex-elementary).

::: function
int acf_exp(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_exp_pi_i(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_log(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sin(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_cos(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sin_cos(acf_t res1, acf_t res2, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sin_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_cos_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sin_cos_pi(acf_t res1, acf_t res2, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sinh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_cosh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sinh_cosh(acf_t res1, acf_t res2, const acf_t x, slong prec, arf_rnd_t rnd)

Set `res` (or the pair `res1`, `res2`) to the exponential, $e^{\pi i x}$,
the principal logarithm, the sine and cosine (of $x$ and of $\pi x$) and
the hyperbolic sine and cosine of `x`. A nan or infinite part gives nan
parts, as do the arguments beyond the ranges of the real functions (the
exponential beyond $2^{\max(128, 2 \mathit{prec})}$ in the real part,
the trigonometric functions beyond $2^{\max(65536, 4 \mathit{prec})}$);
`log(0)` is $-\infty$. Real and purely imaginary arguments are computed
by the real functions (with $\log(-x) = \log x + \pi i$, $\log(ix) =
\log|x| \pm \pi i/2$).
:::

::: function
int acf_pow(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)

Sets `res` to $x^y = \exp(y \log x)$ with the principal logarithm.
$x^0 = 1$ (including $0^0$); $0^y = 0$ for $\operatorname{Re}(y) > 0$
and nan otherwise. Real powers of positive reals are computed by
`arf_pow`; the powers of real and purely imaginary numbers with real
exponents (or with $|x| = 1$) from $|x|^{\operatorname{Re} y}$ and the
sine and cosine of $\pi k \operatorname{Re}(y) / 2$, so that their
exact values (such as $(-2)^3$ and $i^{1000}$) are recognised; see
[](algorithms/special-functions/complex-elementary#acf) for the other exact cases.
:::

::: function
int acf_tan(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_cot(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sec(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_csc(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_tan_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_cot_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sec_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_csc_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_tanh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_coth(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sech(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_csch(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)

Set `res` to the tangent, cotangent, secant and cosecant of `x`
(respectively of $\pi x$) and to the hyperbolic tangent, cotangent,
secant and cosecant of `x`. Real arguments are computed by the real
functions (whose results at the poles are kept, with a zero imaginary
part); a nan or infinite part gives nan parts.

:::

::: function
int acf_expm1(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_log1p(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sinc(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_sinc_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_asin(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_acos(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_atan(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_acot(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_asec(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_acsc(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_asin_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_acos_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_atan_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_asinh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_acosh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_atanh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_acoth(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_asech(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_acsch(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)

Set `res` to $e^x - 1$, $\log(1 + x)$, $\sin(x)/x$, $\sin(\pi x)/(\pi
x)$, and the principal inverse trigonometric and hyperbolic functions
of `x` (and the inverse trigonometric functions divided by $\pi$), with
the branch cuts of [](acb#inverse-trig) and [](acb#inverse-hyperbolic); the reciprocal functions are
$\operatorname{acot}(z) = \operatorname{atan}(1/z)$ and so on. At the
poles and the singular points (such as $\operatorname{atanh}(1)$ and
$\operatorname{acot}(0)$) the parts are nan.

:::

::: function
int acf_inv(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_div(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)

Set `res` to $1/x$, respectively $x/y$, with each part correctly rounded
(the relaxed modes as described above): $1/(a + bi) = (a -
bi)/(a^2 + b^2)$ and $(a + bi)/(c + di) = ((ac + bd) + (bc - ad) i) /
(c^2 + d^2)$ with the numerators and the denominator formed exactly, so
that each part is one correctly rounded quotient of exact numbers; when
the terms of a sum lie more than $2 \mathit{prec} + 128$ bits apart (so
that they cannot cancel) the sum is rounded to `prec` + 64 bits instead,
the quotient certified from its error bound, and a part that this leaves
undecided decided exactly by the sign of $N - \beta s$ for its exact
numerator $N$, the denominator $s$ and the rounding boundary $\beta$
next to it (a sum of exact terms, whatever their sizes).
A real or purely imaginary divisor is handled by real divisions of the
parts; a zero divisor gives the results of `arf_div` (infinities and
nan).
:::

::: function
int acf_sqrt(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
int acf_rsqrt(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)

Set `res` to the principal square root of `x`, respectively its
reciprocal, with each part correctly rounded (the relaxed modes as
described above). The dyadic roots are computed exactly; otherwise the
complex kernels of [](arf#mpn) (or, for exponents beyond them, the
compositions $q = \sqrt{(|x| + |a|)/2}$, $|b|/(2q)$ and their quotients
by $|x|$ at `prec` + 64 bits), and a part that these leave undecided is
decided exactly from the sign of a polynomial in $a$, $b$ and the
rounding boundary next to it, without retries. Real and purely imaginary arguments are handled by the real
roots; `rsqrt(0)` is $+\infty$.
:::
