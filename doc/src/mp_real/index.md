---
title: mp_real.h
header: mp_real.h
description: Limb-radix multiple-precision real balls and fixed-point numbers, with the elementary functions, constants and hypergeometric series built on them.
category: Floating-point arithmetic
order: 3
---

This module provides multiple-precision real arithmetic in the radix
$B = 2^{\bits}$, with limb granularity throughout and as little overhead
per operation as possible. It is the engine behind the high-precision
constants of `arb` (`arb_const_pi`, `arb_const_euler`, ...), the
elementary functions of `nfloat`, the partition function and the
Newton divisions of `mpn_extras`, and it can be used directly as a
backend for numerical algorithms: binary splitting, AGM iterations,
Newton iterations and evaluation of elementary functions. It works at
two levels.

- `mp_real_t` is a ball: a midpoint with a limb mantissa and a
  radix-$B$ exponent, and a rigorous radius counted in units of the
  mantissa's last limb. The functions `mp_real_*` provide arithmetic,
  roots, the arithmetic-geometric mean, the elementary functions,
  constants and hypergeometric series on balls.
- Fixed-point numbers in limb arrays represent values in $[0, 1)$ (or
  with one or two integral limbs). The functions `_mp_real_*` on them
  provide the kernels of the elementary functions with explicit error
  bounds in ulps, Newton inverses and square roots, and exact floors of
  constants. A leading underscore also marks ball functions with
  unchecked preconditions.

The module is tuned for 64-bit machines; with 32-bit limbs the
generated straight-line kernels are disabled and generic code paths are
used.

The guides below explain the two levels: [](mp_real/balls) describes
the ball type, precision and radius, and [](mp_real/fixed-point) the
fixed-point kernels and their error contract. The algorithms are
described in [](algorithms/floating-point/mp-real-balls),
[](algorithms/floating-point/elementary-functions),
[](algorithms/floating-point/diophantine-reduction),
[](algorithms/floating-point/hypergeometric-series) and
[](algorithms/floating-point/constants).

## Types and conventions {#types}

::: type
mp_real_struct
mp_real_t

A ball

$$
(-1)^{\mathit{negative}}\, (d, \mathit{size})\, B^{\mathit{exp} - \mathit{size}}
\;\pm\; \mathit{err}\, B^{\mathit{exp} - \mathit{size}}
$$

with fields `d` (the mantissa limbs, least significant first), `alloc`,
`size`, `negative`, `exp` (an `slong` exponent in limbs) and `err` (a
`ulong` radius in ulps of the mantissa's last limb). A nonzero mantissa
has a nonzero top limb, so that $B^{\mathit{exp}-1} \le |\mathrm{mid}| <
B^{\mathit{exp}}$; `size` = 0 encodes a zero midpoint, whose radius is
counted in units of $B^{\mathit{exp}}$. An `mp_real_t` is an array of
one `mp_real_struct`, so that it is passed by reference.

::: details
The radius count satisfies $0 \le \mathit{err} < B$, and `err` = 0 means
that the ball is exact. The radius has no exponent of its own; instead
the mantissa length is adjusted to it:

- a radius below one ulp of the current mantissa is installed by padding
  the mantissa with zero limbs down to the radius's scale (so balls such
  as $1 - \varepsilon$ carry their zero limbs explicitly);
- a radius of $B$ ulps or more truncates the mantissa by as many limbs
  as needed, the dropped limbs being rounded into one more unit at the
  new last limb.

The mantissa length therefore tracks the number of accurate limbs, with
at most one limb of noise. Exact values carry no low zero limbs, so that
small exact integers stay one limb long. See
[](algorithms/floating-point/mp-real-balls) for the consequences.
:::
:::

::: type
mp_real_ptr
mp_real_srcptr

Pointers to `mp_real_struct` and `const mp_real_struct`, used for
vectors of balls.
:::

A *fixed-point number* $(x, n)$ is an array of $n$ limbs
$x_0, \ldots, x_{n-1}$ representing the fraction
$\sum_i x_i B^{i-n} \in [0, 1)$, with unit in the last place (ulp)
$B^{-n}$. An output written to $n + 1$ limbs has an additional integral
(units) limb at index $n$, and one written to $n + 2$ limbs two integral
limbs. The fixed-point functions with an ulp-accurate contract take an
argument `ulong * err` immediately after their outputs; they write there
a bound for $|\mathrm{res} - f(x)|$ in ulps of the output ($B^{-n}$),
valid for each output of a function with two outputs, and `err` may be
`NULL`. Each function documents the largest value it can write, so that
callers can size their guard limbs. Unless stated otherwise the bounds
are two-sided. The guide [](mp_real/fixed-point) gives more detail.

All precisions $n$ are numbers of limbs, except that the elementary
functions on balls take a precision in bits (`_bits` suffix).

::: summary algorithms/floating-point/mp-real-balls :::

## Exponent safety and working precision {#precision}

The ball arithmetic is rigorous except that exponents are not checked
for overflow; a caller that can produce huge exponents checks its
operands against a safe range instead.

::: macro
MP_REAL_EXP_MAX

The largest safe limb exponent, `WORD_MAX / (4 * FLINT_BITS)`. Since
several functions form the bit exponent `FLINT_BITS * exp`, this keeps
bit exponents below `WORD_MAX / 4`, which is $2^{61}$ on 64-bit
machines.
:::

::: function
int mp_real_exp_is_safe(slong e)
int mp_real_is_safe(const mp_real_t x)

Return whether $|e| \le$ `MP_REAL_EXP_MAX`, respectively whether the
exponent of $x$ is safe.

::: details
An operation whose ball operands are safe, and whose bit-exponent
arguments $e$ (of `mp_real_mul_2exp_si`, `mp_real_add_error_2exp_si`,
`mp_real_add_rel_error_2exp_si` and `_mp_real_set_mpn_2exp`) satisfy
$|e| \le$ `FLINT_BITS * MP_REAL_EXP_MAX`, computes without overflow. Its
result may lie outside the safe range, so a chain of operations whose
exponents can grow must check again.
:::
:::

::: macro
MP_REAL_PREC_GUARD

The number of guard limbs added by `mp_real_prec_limbs` and
`mp_real_prec_bits`, currently 1.
:::

::: function
slong mp_real_prec_limbs(slong prec_limbs)
slong mp_real_prec_bits(slong prec_bits)

Return the working precision $n$ in limbs for a target relative
accuracy $\varepsilon = B^{-\mathit{prec\_limbs}}$ resp.
$2^{-\mathit{prec\_bits}}$: `prec_limbs` resp.
$\lceil \mathit{prec\_bits} / \bits \rceil$ plus `MP_REAL_PREC_GUARD`.

::: details
At this precision the basic operations `mp_real_add`, `mp_real_sub`,
`mp_real_mul`, `mp_real_div`, `mp_real_sqrt`, `mp_real_rsqrt`,
`mp_real_mul_ui`, `mp_real_div_ui` and `mp_real_addmul_ui` return a
relative radius of at most $4\varepsilon$ on exact operands and at most
$16\varepsilon$ on operands of relative radius at most $\varepsilon$
(for sums of operands of like signs and differences without
cancellation). Without the guard limb a product of exact operands could
lose a whole limb.
:::
:::

## Memory management and assignment {#memory}

::: function
void mp_real_init(mp_real_t x)
void mp_real_clear(mp_real_t x)

Initialise $x$ to exact zero, respectively free its memory.
:::

::: function
void mp_real_fit_length(mp_real_t x, slong len)
void _mp_real_grow(mp_real_t x, slong len)

Ensure that $x$ has room for `len` mantissa limbs. The inline function
calls `_mp_real_grow` when it has to reallocate; the contents of the
mantissa are preserved.
:::

::: function
void mp_real_swap(mp_real_t x, mp_real_t y)

Swap $x$ and $y$ efficiently.
:::

::: function
void mp_real_zero(mp_real_t x)
void mp_real_set(mp_real_t res, const mp_real_t x)
void mp_real_set_ui(mp_real_t x, ulong c)
void mp_real_set_si(mp_real_t x, slong c)
void mp_real_set_fmpz(mp_real_t x, const fmpz_t f)

Set the ball to exact zero, to a copy of $x$, or to the exact integer
$c$ resp. $f$.

::: details
`mp_real_set_fmpz` copies the limbs of $f$ above its low zero limbs
directly into the mantissa, without temporaries.
:::
:::

::: function
void mp_real_set_trunc(mp_real_t res, const mp_real_t x, slong n)

Set `res` to $x$ rounded to its top $n \ge 1$ limbs (a copy if $x$ has
at most $n$ limbs). `res` may alias $x$.

::: details
The dropped tail and the old radius are each below one ulp of the new
last limb, so the new radius is $[\mathrm{tail} \ne 0] +
[\mathrm{radius} \ne 0]$ ulps. The result is exact if nothing nonzero is
dropped from an exact input.
:::
:::

::: function
void _mp_real_set_mpn_2exp(mp_real_t x, nn_srcptr p, slong len, slong e)

Set $x$ exactly to $(p, \mathit{len})\, 2^e$ for the unsigned integer
$(p, \mathit{len})$, which may have zero top limbs. The bit part of $e$
costs one limb shift.
:::

::: function
void mp_real_neg(mp_real_t res, const mp_real_t x)

Set `res` to $-x$; in place this flips the sign bit.
:::

::: function
void mp_real_mul_2exp_si(mp_real_t res, const mp_real_t x, slong e)

Set `res` to $x \cdot 2^e$ exactly (an exponent change and at most one
limb shift of the mantissa).
:::

## Predicates, magnitudes and the radius {#radius}

::: function
int mp_real_is_zero(const mp_real_t x)

Return whether $x$ is exactly zero (zero midpoint and zero radius).
:::

::: function
slong mp_real_abs_bound_lt_2exp_si(const mp_real_t x)

Return an $e$ with $|t| < 2^e$ for every $t$ in the ball $x$.
:::

::: function
slong mp_real_rel_radius_lt_2exp_si(const mp_real_t x)

Return an $e$ such that the radius of $x$ is below $2^e$ times the
absolute value of its midpoint: `-WORD_MAX / 2` for an exact $x$ and
`WORD_MAX / 2` for a zero midpoint with a nonzero radius.
:::

::: function
void mp_real_add_error_2exp_si(mp_real_t x, slong e)
void mp_real_add_rel_error_2exp_si(mp_real_t x, slong e)

Add $2^e$ to the radius of $x$, respectively enlarge $x$ to contain
$x(1 + t)$ for every $|t| \le 2^e$.
:::

::: function
void _mp_real_add_error_ulps(mp_real_t x, double v)
void _mp_real_add_error_ulps_at(mp_real_t x, double v, slong anchor)

Add $v$ ulps of the mantissa's last limb, respectively $v
B^{\mathit{anchor}}$, to the radius of $x$, for $v \ge 0$.

::: details
The first form depends on the current normalisation of the mantissa:
an exact value has no low zero limbs, so its ulp can be much coarser
than the scale the caller has in mind. The second form states the scale
explicitly and should be preferred.
:::
:::

## Arithmetic {#arithmetic}

These functions round their outputs to about $n$ limbs, or fewer when an
operand is accurate to fewer, and add all rounding errors to the radius.
Outputs may alias inputs. The guide [](mp_real/balls) explains how to
choose $n$.

::: function
void mp_real_add(mp_real_t res, const mp_real_t a, const mp_real_t b, slong n)
void mp_real_sub(mp_real_t res, const mp_real_t a, const mp_real_t b, slong n)

Set `res` to $a + b$ resp. $a - b$.

::: details
The sum is formed in the destination without a temporary. A sum taken
in place, `x := x + y` with $y$ inside the window of $x$, costs the limb
addition of $y$ plus a constant, independently of the length of $x$. A
$y$ reaching below the last limb of an inexact $x$ is cut there, its
tail entering the radius. Only a $y$ below an exact $x$ (the contiguous
mantissa must move) or an $x$ longer than the precision asks for (a
truncation) cost a pass over $x$.
:::
:::

::: function
void mp_real_addmul_ui(mp_real_t res, const mp_real_t x, const mp_real_t y, ulong c, slong n)
void mp_real_submul_ui(mp_real_t res, const mp_real_t x, const mp_real_t y, ulong c, slong n)

Set `res` to $x + yc$ resp. $x - yc$ for a word $c$.

::: details
In place (`res` aliasing $x$) and with $yc$ inside the window of $x$
(and, for a difference, unable to flip its sign), this is one
`mpn_addmul_1` resp. `mpn_submul_1` over the limbs of $y$, independent
of the length of $x$; $y$ may reach below an inexact $x$, the tail
entering the radius. These are the primitives for summing many terms
with word coefficients into one accumulator, as in rectangular
splitting. Other shapes form $yc$ and add it.
:::
:::

::: function
void mp_real_mul(mp_real_t res, const mp_real_t a, const mp_real_t b, slong n)
void mp_real_mul_ui(mp_real_t res, const mp_real_t a, ulong c, slong n)

Set `res` to $ab$ resp. $ac$.

::: details
The operands are first truncated to $n + 2$ limbs, the dropped tails
entering the radius through the other operand's magnitude, so that a
long exact operand costs no more than the precision asks for. Short
products are formed in full and longer ones as windowed middle products
(`flint_mpn_mulmid`); a square goes through `flint_mpn_sqr`.
:::
:::

::: function
void mp_real_mul_complex(mp_real_t rr, mp_real_t ri, const mp_real_t ar, const mp_real_t ai, const mp_real_t br, const mp_real_t bi, slong n)

Set $r_r + i r_i$ to $(a_r + i a_i)(b_r + i b_i)$, to about $n$ limbs
relative to the larger part of the result. The outputs may alias the
inputs but not each other.

::: details
Both parts of the result, and both parts of each operand, are kept down
to the same last limb: the precision of a complex number is that of its
modulus, and an accumulation of $\prod_k \exp(i x_k)$ needs its small
imaginary part to the absolute precision of the frame. The product is
one call to `flint_mpn_mul_complex` (its high half,
`flint_mpn_mulhigh_n_complex`, for long operands), which shares the
transforms of the four parts instead of making four products.
:::
:::

::: function
void _mp_real_submul_bounded(mp_real_t res, const mp_real_t a, const mp_real_t b, const mp_real_t c, slong E, slong n)

Set `res` to $a - bc$, given that $|a - bc| < B^E$. Only the window of
the product that reaches the result is computed, and the result is wrong
if the bound does not hold.

::: details
This is the residual of a Newton or Taylor step, whose bound comes from
the radius of the previous approximation. The sign of $a - bc$ is read
off the top limb of the product window, so the bound must be rigorous.
:::
:::

::: function
void mp_real_div(mp_real_t res, const mp_real_t a, const mp_real_t b, slong n)
void mp_real_div_ui(mp_real_t res, const mp_real_t a, ulong c, slong n)

Set `res` to $a / b$ resp. $a / c$. The divisor ball must have a
relative radius below $2^{-30}$; `mp_real_div` throws for a wider or
zero divisor, and `mp_real_div_ui` for $c = 0$.

::: details
A divisor ball with a relative radius of $2^{-30}$ or more is treated as
a usage error, since the quotient's mantissa would be noise; use
`mp_real_get_arb` and `arb_div` for such balls. `mp_real_div_ui`
divides by `mpn_divrem_1` over appended zero limbs, keeping the
quotient's own precision (an exact dividend gets a quotient of $n + 2$
limbs), so that a chain of divisions by a common denominator loses no
accuracy relative to the quotients.
:::
:::

## Roots and the arithmetic-geometric mean {#roots}

::: function
void mp_real_sqrt(mp_real_t res, const mp_real_t x, slong n)
void mp_real_rsqrt(mp_real_t res, const mp_real_t x, slong n)

Set `res` to $\sqrt{x}$ resp. $1/\sqrt{x}$ for a positive ball $x$ of
relative radius below $2^{-30}$. A wider ball, or a ball with a
negative or zero midpoint, throws, except that `mp_real_sqrt` of an
exact zero gives zero. `res` may alias $x$.

::: details
The exponent is aligned to an even number of limbs, the root computed by
`_mp_real_sqrt_newton` resp. `_mp_real_rsqrt_newton` (the square root
by `flint_mpn_sqrtrem` below a cutoff of a few hundred limbs), and the
radius of $x$ propagated through the derivative.
:::
:::

::: function
void mp_real_rsqrt_ui(mp_real_t res, ulong c, slong n)

Set `res` to $1/\sqrt{c}$ for a word $c \ge 1$, by
`_mp_real_rsqrt_ui_newton`; throws for $c = 0$.
:::

::: macro
MP_REAL_ROOT_K_MAX

The bound on $k$ for `mp_real_root_ui` and `mp_real_rroot_ui`: $2^{40}$ on 64-bit
machines and $2^{27}$ on 32-bit machines.
:::

::: function
void mp_real_root_ui(mp_real_t res, const mp_real_t x, ulong k, slong n)
void mp_real_rroot_ui(mp_real_t res, const mp_real_t x, ulong k, slong n)
void _mp_real_root_ui_order(mp_real_t res, const mp_real_t x, ulong k, slong n, int r, int recip)

Set `res` to $x^{1/k}$ resp. $x^{-1/k}$ for $x > 0$ and $1 \le k <$
`MP_REAL_ROOT_K_MAX`. They throw for $k = 0$, for $x$ not positive (the
third function for $k < 2$ as well). The third function takes the
order $r$ of the iteration ($2 \le r \le 15$, or 0 for the default) and
`recip` = 1 for the reciprocal root.

::: details
For $k = 1$ the result is $x$ resp. $1/x$, and for $k = 2$ these call
`mp_real_sqrt` and `mp_real_rsqrt`. Otherwise one step of order $r$ of a
binomial-series iteration is taken from an approximation $z$ of
$x^{-1/k}$ accurate to about $1/r$ of the precision, obtained
recursively (the recursion starts from a double):

$$
x^{-1/k} = z (1 - u)^{-1/k}, \qquad
x^{1/k} = x z^{k-1} (1 - u)^{-(k-1)/k}, \qquad u = 1 - x z^k,
$$

with the series $(1 - u)^{-b} = \sum_j c_j u^j$ truncated after
$u^{r-1}$ and its tail, below $2 c_r |u|^r$ for $|u| \le 1/2$, added to
the radius. The series is summed over a common denominator of the
$c_j$, so that it costs one division by an integer; the coefficient
tables are cached per thread. The default order is 3 for $k < 8$ and 4
otherwise. Only the size of the final $u$ is computed rigorously, so
the radius of the recursive approximation is discarded.
:::
:::

::: function
void mp_real_agm(mp_real_t res, const mp_real_t x, const mp_real_t y, slong n)
void _mp_real_agm_order(mp_real_t res, const mp_real_t x, const mp_real_t y, slong n, int m)

Set `res` to the arithmetic-geometric mean $\operatorname{agm}(x, y)$
of $x, y \ge 0$ to about $n$ limbs; throws for a negative midpoint. The
second function takes the order $m$ of the final series (at most 16, 10
on 32-bit machines; 0 for the default).

::: details
The iteration $a' = (a + b)/2$, $b' = \sqrt{ab}$ runs in ball arithmetic
until $z = (a - b)/(a + b)$ satisfies $z^2 < 2^{-p/m}$ for $p$ bits, and
is finished by

$$
\operatorname{agm}(a, b) = \frac{(a + b)/2}{{}_2F_1(\tfrac12, \tfrac12; 1; z^2)},
\qquad
\frac{1}{{}_2F_1(\tfrac12, \tfrac12; 1; t)} = 1 - \sum_{j \ge 1} c_j t^j.
$$

The $c_j$ are positive dyadic rationals summing to 1 (the reciprocal
vanishes at $t = 1$; positivity follows from the log-convexity of the
coefficients of ${}_2F_1$ by Kaluza's theorem), so the tail after
$t^{m-1}$ is below $t^m$. The default $m$ is 8 below 1024 limbs and 12
above. When $a - b$ no longer determines its own size, or $a + b$ or
$ab$ is too wide for the division and the square root, the result is
the enclosure $(a + b)/2 \pm |a - b|/2$, valid because
$b \le \operatorname{agm}(a, b) \le a$. The working set is four
numbers of $n$ limbs.
:::
:::

## Conversions and output {#conversions}

::: function
void mp_real_get_arb(arb_t res, const mp_real_t x)

Set `res` to the ball $x$ without loss: the midpoint exactly and the
radius rounded up.
:::

::: function
int mp_real_unique_integer(mp_real_t res, const mp_real_t x)

If the ball $x$ contains exactly one integer, set `res` to that integer
(exactly) and return 1; otherwise return 0 and leave `res` unchanged.
`res` may alias $x$, in which case no reallocation takes place.

::: details
This is `arb_get_unique_fmpz` with the result kept as an `mp_real_t`,
free of the size limits of `fmpz_t`. With a fractional part of $f$ limbs and
a radius below one ulp, the candidates follow from the carry of the
fraction plus the radius and the borrow of the fraction minus the
radius, decided by scanning the fraction from its top limb until a limb
is neither all zeros nor all ones; the integer part is then moved down
(and incremented if needed).
:::
:::

::: function
void mp_real_set_dfloat(mp_real_t res, const double * x, slong n, double rad)
void mp_real_get_dfloat(double * res, double * rad, slong n, const mp_real_t x)

Convert from and to the double expansions of `dfloat`: arrays of
$n \le 8$ doubles (the midpoint components) and a double radius.
`mp_real_set_dfloat` sets `res` to the ball $\sum_i x_i \pm \mathit{rad}$
and throws for nonfinite components or radius. `mp_real_get_dfloat`
rounds $x$ greedily to $n$ doubles and sets `rad` (if not `NULL`) to a
bound for the error and the radius of $x$.

::: details
The sum in `mp_real_set_dfloat` is exact, accumulated in a two's
complement frame spanning the components' exponents, and the radius is
installed as an ulp count, rounded up. In `mp_real_get_dfloat` each
component is the double nearest to the remainder, so that
$|x_{i+1}| \le \operatorname{ulp}(x_i)/2$; components that would be
subnormal are not produced (their value goes to the radius), and a value
beyond the double range gives dfloat's whole line, zero components with
an infinite radius. Both conversions avoid `arb_t`.
:::
:::

::: function
void mp_real_print(const mp_real_t x)

Print the fields `size`, `exp` and `err` of $x$ and its value to about
twenty digits.
:::

::: function
void _mp_real_get_fixed(nn_ptr y, ulong * err, const mp_real_t x, slong n)

Write a ball known to lie in $[0, 1)$ as the fixed-point number
$(y, n)$, truncating, and set `err` (if not `NULL`) to a bound for the
error in ulps $B^{-n}$, rounded up; `UWORD_MAX` means that no bound
below $B - 1$ is known. A negative midpoint within the radius is clamped
to zero and absorbed into the bound.
:::

::: function
int _mp_real_get_fixed_floor(nn_ptr y, slong n, const mp_real_t x)

Write $y = \lfloor x B^n \rfloor$ for a ball known to lie in $[0, 1)$
and return 1 if the radius determines that floor; otherwise return 0, so
that the caller can retry at a higher precision.
:::

## Constants {#constants}

The constants are computed by binary splitting of hypergeometric
series ([](algorithms/floating-point/constants)) and can be cached
process-wide. [](mp_real/constants) describes the functions in
detail.

::: function
void mp_real_const_pi4(mp_real_t res, slong n, int cache)
void mp_real_const_log2(mp_real_t res, slong n, int cache)
void mp_real_const_euler(mp_real_t res, slong n, int cache)
void mp_real_const_e(mp_real_t res, slong n, int cache)
void mp_real_const_log10(mp_real_t res, slong n, int cache)
void mp_real_const_catalan(mp_real_t res, slong n, int cache)
void mp_real_const_zeta3(mp_real_t res, slong n, int cache)
void mp_real_const_zeta5(mp_real_t res, slong n, int cache)
void mp_real_const_gamma_1_3(mp_real_t res, slong n, int cache)
void mp_real_const_gamma_1_4(mp_real_t res, slong n, int cache)
void mp_real_const_2_div_pi(mp_real_t res, slong n, int cache)
void mp_real_const_inv_log2_frac(mp_real_t res, slong n, int cache)

Set `res` to a ball containing $c = \pi/4$, $\log 2$, Euler's constant
$\gamma$, $e$, $\log 10$, Catalan's constant $G$, $\zeta(3)$,
$\zeta(5)$, $\Gamma(1/3)$, $\Gamma(1/4)$, $2/\pi$ resp.
$1/\log 2 - 1$, accurate to about $n$ limbs. With `cache` = 0 the value
is computed from scratch; with `cache` = 1 it is read from (and if
necessary added to) a process-wide cache of floors, and the ball is the
cached floor with a radius of one ulp.

::: details
The cache holds, per constant, $\lfloor c B^N \rfloor$ for the largest
$N$ computed so far by any thread. Floors nest (the top limbs of
$\lfloor c B^N \rfloor$ are $\lfloor c B^n \rfloor$ for $n \le N$), so a
shorter request is a copy, and a longer one recomputes the entry at
$\max(n + 5, 1.5 N)$ limbs, so that a sequence of growing requests costs
a constant factor more than the last one. $\pi/4$, $\log 2$, $2/\pi$ and
$1/\log 2 - 1$ are read from static tables up to 64 limbs. The last two
serve the argument reductions modulo $\pi/2$ and $\log 2$.
:::
:::

::: function
void _mp_real_const_pi4(nn_ptr res, ulong * err, slong n, int cache)
void _mp_real_const_log2(nn_ptr res, ulong * err, slong n, int cache)
void _mp_real_const_euler(nn_ptr res, ulong * err, slong n, int cache)
void _mp_real_const_e(nn_ptr res, ulong * err, slong n, int cache)
void _mp_real_const_log10(nn_ptr res, ulong * err, slong n, int cache)
void _mp_real_const_catalan(nn_ptr res, ulong * err, slong n, int cache)
void _mp_real_const_zeta3(nn_ptr res, ulong * err, slong n, int cache)
void _mp_real_const_zeta5(nn_ptr res, ulong * err, slong n, int cache)
void _mp_real_const_gamma_1_3(nn_ptr res, ulong * err, slong n, int cache)
void _mp_real_const_gamma_1_4(nn_ptr res, ulong * err, slong n, int cache)
void _mp_real_const_2_div_pi(nn_ptr res, ulong * err, slong n, int cache)
void _mp_real_const_inv_log2_frac(nn_ptr res, ulong * err, slong n, int cache)

Set `res` to exactly $\lfloor c B^n \rfloor$ for the constants $c$ of
`mp_real_const_pi4` and its siblings,
and `err` (if not `NULL`) to 1. The output has $n$ limbs for the
constants below 1 ($\pi/4$, $\log 2$, $\gamma$, $G$, $2/\pi$,
$1/\log 2 - 1$) and $n + 1$ limbs, with a units limb, for the others.
`cache` is as for the ball versions.

::: details
The floor is read off a ball computed with a few guard limbs, which is
recomputed with more guard limbs until its radius determines the floor.
The value 1 written to `err` bounds $c - \lfloor c B^n \rfloor B^{-n}$
in ulps; the error is one-sided.
:::
:::

::: function
void _mp_real_const_clear_cache(void)

Free the shared caches of the constants. No other thread may be using
them. `flint_cleanup_master` also frees them.
:::

::: type
mp_real_const_func

The type `void (*)(mp_real_t, slong, int)` of the ball constant
functions (`mp_real_const_pi4` and so on).
:::

::: function
void _mp_real_const_arb(arb_t res, mp_real_const_func f, slong prec)

Set `res` to the constant computed by `f` with `cache` = 0, as an `arb`
ball of `prec` bits. This is the entry point of the `arb` constants
(`arb_const_pi`, `arb_const_log2`, `arb_const_euler`, ...), which keep
their own caches.
:::

::: function
void _mp_real_const_euler_tune(mp_real_t res, slong n, int set)

Euler's constant with a forced choice of the logarithms used in the
Brent–McMillan formula: `set` = 0, 1, 2, 3 selects the primes $\{2\}$,
$\{2, 3\}$, $\{2, 3, 5\}$, $\{2, 3, 5, 7\}$, and $-1$ the automatic
choice. For tuning.
:::

## Elementary functions {#elementary}

These functions evaluate the elementary functions on balls to a
relative accuracy of about $2^{-\mathit{prec}}$, or less as the radius
of the argument allows. The midpoint is evaluated exactly as given, the
working precision is chosen internally, and the outputs may alias the
argument. [](mp_real/elementary) explains the common
behaviour, and [](algorithms/floating-point/elementary-functions) the
algorithms.

::: function
void mp_real_sin_cos_bits(mp_real_t res1, mp_real_t res2, const mp_real_t x, slong prec)

Set `res1` and `res2` to balls containing $\sin y$ and $\cos y$ for
every $y$ in $x$. Either output may be `NULL`. The sine keeps its
relative accuracy for small $|x|$; near the other zeros the loss of
significance shows in the radius. Arguments with $|x| \ge
2^{\max(65536, 4\,\mathit{prec})}$ give $[0 \pm 1]$.

::: details
If $x$ is inexact, `prec` is first lowered to a few bits beyond the
accuracy its radius allows (relative for $|x| < 1$, absolute
otherwise), and the radius is added to both outputs at the end (as
$|\sin'|, |\cos'| \le 1$). A ball inside $(-1, 1)$ whose radius exceeds
a quarter of its midpoint gives $[0 \pm 2^e]$ and $1 \pm 2^{2e-1}$ for
$|x| < 2^e$; a wider ball elsewhere gives $[0 \pm 1]$.

For a midpoint with $z$ leading zero bits ($|m| < 2^{-z}$): if
$2z \ge \mathit{prec} + 3$ the results are $m$ and 1, and if $4z \ge
\mathit{prec} + 3$ they are $m - m^3/6$ and $1 - m^2/2$, with the Taylor
remainders in the radius. Otherwise, for $|m| < 1$, the fixed-point
kernel is called on $|m|$ at $\mathit{prec} + z$ bits plus guard bits
(or a series for the reduced argument when $z$ is large). For $|m| \ge
1$ the argument is reduced modulo $\pi/2$, with the quotient chosen so
that the remainder is nonnegative, folded to $[0, \pi/4]$, and the
outputs are permuted and negated according to the quotient modulo 4.
The reduction adds at most 3 ulps to the kernel's bound and is not
repeated at a higher precision near zeros.
:::
:::

::: function
int mp_real_tan_bits(mp_real_t res, const mp_real_t x, slong prec)

Set `res` to a ball containing $\tan y$ for every $y$ in $x$ and return
1, or return 0 with `res` set to zero when the ball may contain a pole
(or $|x| \ge 2^{\max(65536, 4\,\mathit{prec})}$). The result keeps its
relative accuracy near the zeros and poles of the tangent.

::: details
After the reduction $|m| = a\pi/2 + sv$, $v \in [0, \pi/4]$, the result
is $s \tan v$ for even $a$ and $-s/\tan v$ for odd $a$ (one Newton
reciprocal). Near a zero or a pole $v$ is small, so the kernel runs with
the leading zero bits of $v$ as extra precision, and the reduction is
repeated at that precision when the first pass had fewer limbs. An exact
$m$ is never a multiple of $\pi/2$, and by the irrationality measure of
$\pi$ a $b$-bit $m$ lies at least about $2^{-7.2 b}$ away from one; the
repetition is capped at $8b + 1024$ extra bits (status 0 beyond). For
an inexact $x$ no repetition is made once $v$ falls below the radius.

The radius $r$ of $x$ enters through a pole test: with $T$ enclosing
$\tan m$ and $U = r(1 + |T|) \ge r/|\cos m|$, $U < 1$ implies that the
ball contains no pole, and then
$|\tan y - \tan m| \le r(1 + T^2)/(1 - U)^2$.
:::
:::

::: function
void mp_real_sin_cos_pi_bits(mp_real_t res1, mp_real_t res2, const mp_real_t x, slong prec)
int mp_real_tan_pi_bits(mp_real_t res, const mp_real_t x, slong prec)

The same as `mp_real_sin_cos_bits` and `mp_real_tan_bits` for
$\sin(\pi y)$, $\cos(\pi y)$ and $\tan(\pi y)$. The argument reduction
is exact, so the outputs keep their relative accuracy near every zero
(and the tangent near every pole) at any magnitude of $x$, and the
values at multiples of $1/2$ (and of the tangent at odd multiples of
$1/4$) are exact.

::: details
The midpoint is dyadic, and its units bit and fraction bits give
$|m| = a/2 + st$ with $a \in \{0, 1, 2, 3\}$, $s = \pm 1$ and
$t \in [0, 1/4]$ exactly. Then $v = \pi t$ enters the formulas of the
reduction modulo $\pi/2$; $t = 0$ gives $0$ or $\pm 1$ (and a pole of
the tangent for odd $a$, status 0), and $t = 1/4$ the exact
$\tan(\pi/4) = 1$. Otherwise $t$ has $z$ leading zero bits and the
kernel runs at $z$ extra bits. The radius $r$ of $x$ is a radius $\pi r$
of the argument; $\pi r \ge 1/4$ gives $[0 \pm 1]$ for the sine and
cosine.
:::
:::

::: function
void mp_real_sin_cos_pi_ui_div_ui(mp_real_t res1, mp_real_t res2, ulong p, ulong q, slong n)
int mp_real_tan_pi_ui_div_ui(mp_real_t res, ulong p, ulong q, slong n)

Set `res1` and `res2` to balls containing $\sin(\pi p/q)$ and
$\cos(\pi p/q)$, respectively `res` to $\tan(\pi p/q)$, for words
$p \ge 0$ and $q \ge 1$, to a relative accuracy of about $B^{-n}$. Either
output of the first function may be `NULL`; the two must be distinct.
The tangent returns 0 with `res` set to zero at a pole, else 1. The
rational values $0, \pm 1, \pm 1/2$ are exact. Both throw for $q = 0$.

::: details
The fraction is reduced to lowest terms and folded by the symmetries to
$\pi t/q$ with $0 \le 2t \le q$. The nonzero values of $\sin(\pi t/q)$
and $\cos(\pi t/q)$ are at least $1/q$, so
`_mp_real_sin_cos_pi_ui_div_ui` at
$n + \lceil(\operatorname{bits}(q) + 16)/\bits\rceil$ limbs gives the
relative accuracy. The tangent uses the closed forms $\sqrt 3$,
$1/\sqrt 3$, $\sqrt 2 \mp 1$, $2 \mp \sqrt 3$ for $q = 3, 6, 8, 12$ and
$\sqrt{5 \mp 2\sqrt 5}$, $\sqrt{1 \mp 2/\sqrt 5}$ for $q = 5, 10$ (from
12 limbs), the Chebyshev sine and cosine and one division for the
denominators where they win, and the tangent kernel otherwise.
:::
:::

::: function
void mp_real_exp_bits(mp_real_t res, const mp_real_t x, slong prec)
int mp_real_log_bits(mp_real_t res, const mp_real_t x, slong prec)
void mp_real_atan_bits(mp_real_t res, const mp_real_t x, slong prec)

Set `res` to a ball containing $\exp(y)$, $\log(y)$ resp.
$\operatorname{atan}(y)$ for every $y$ in $x$. The results keep their
relative accuracy near $\exp(y) = 1$, $\log(y) = 0$ and
$\operatorname{atan}(y) = 0$. `mp_real_exp_bits` throws for
$|y| \ge 2^{\bits - 5}$, where the exponent of the result would leave
the safe range. `mp_real_log_bits` returns 0 with `res` set to zero
unless the ball is strictly positive, and 1 otherwise.

::: details
The radius $\rho$ of $x$ enters as the factor $[1 \pm (e^\rho - 1)]$
(exp), as $\rho/(m - \rho)$ (log, by the mean value theorem), and as
$\rho$, or $\rho/(|m| - \rho)^2$ when the ball stays beyond 1 in
absolute value (atan). An exp ball too wide for two correct bits gives
$[0 \pm e^{m + \rho}]$.

- exp: $|m| = q \log 2 + t$, $t \in [0, \log 2)$, with $q$ from two
  limbs of $1/\log 2$ and one multiply-subtract by $\log 2$; for $m < 0$
  the quotient is rounded up and the remainder complemented, so that
  $\exp(m) = 2^{-q}\exp(t)$ directly. Small $|m|$ use $1 + m$,
  $1 + m + m^2/2$, ... or a series for the reduced argument.
- log: with $m = 2^E u$, $u \in [1/2, 1)$, $\log m = (E - 1)\log 2 +
  \operatorname{log1p}(2u - 1)$ up to 600 limbs and $E \log 2 - (-\log u)$
  by `_mp_real_neglog_newton` above; near 1, $2\operatorname{atanh}(s/(2
  \pm s))$ for $s = |m - 1|$, or the first terms of $\log(1 + D)$.
- atan: $\operatorname{atan}(m) = \pi/4 + \operatorname{atan}(s/(2+s))$
  for $1 < m < 2$, $s = m - 1$, and $\pi/2 - \operatorname{atan}(1/m)$
  for $m \ge 2$; $\operatorname{atan}(1) = \pi/4$.

At few limbs the per-size kernels are called directly on the argument's
limbs in place.
:::
:::

::: function
int mp_real_log1p_bits(mp_real_t res, const mp_real_t x, slong prec)
void mp_real_expm1_bits(mp_real_t res, const mp_real_t x, slong prec)
void mp_real_sinh_cosh_bits(mp_real_t rs, mp_real_t rc, const mp_real_t x, slong prec)
void mp_real_tanh_bits(mp_real_t res, const mp_real_t x, slong prec)
int mp_real_asin_bits(mp_real_t res, const mp_real_t x, slong prec)
int mp_real_acos_bits(mp_real_t res, const mp_real_t x, slong prec)
void mp_real_asinh_bits(mp_real_t res, const mp_real_t x, slong prec)
int mp_real_acosh_bits(mp_real_t res, const mp_real_t x, slong prec)
int mp_real_atanh_bits(mp_real_t res, const mp_real_t x, slong prec)

Set `res` (`rs` to $\sinh$ and `rc` to $\cosh$, either may be `NULL`)
to a ball containing $\log(1 + y)$, $e^y - 1$, $\sinh y$, $\cosh y$,
$\tanh y$, $\operatorname{asin} y$, $\operatorname{acos} y$,
$\operatorname{asinh} y$, $\operatorname{acosh} y$ resp.
$\operatorname{atanh} y$ for every $y$ in $x$. The results keep their
relative accuracy at the zeros of the functions, and tiny and huge
arguments are handled without squaring them. The functions returning
`int` return 0 with `res` set to zero when the ball is not certified to
lie in the domain, and 1 otherwise.

::: details
The domains are: $1 + y > 0$ for log1p, $|y| < 1$ for atanh, and for
asin, acos and acosh $|y| \le 1$ resp. $y \ge 1$ for an exact argument
and strictly inside for an inexact one. `mp_real_expm1_bits` throws like
`mp_real_exp_bits` for $y \ge 2^{\bits - 5}$, and
`mp_real_sinh_cosh_bits` for $|y| \ge 2^{\bits - 5}$; expm1 of a large
negative argument gives $-1$ with a tiny radius.

The functions are compositions of `mp_real_exp_bits`,
`mp_real_log_bits` and `mp_real_atan_bits` by formulas without
cancellation:
$\operatorname{log1p}(y) = \log(1 + y)$ with $1 + y$ formed exactly (or
to $\mathit{prec} + z$ bits for $|y| < 2^{-z}$);
$\operatorname{expm1}$ from $\exp$ at $\mathit{prec} + z$ bits;
$\sinh$ and $\cosh$ from $\operatorname{expm1}(\pm m)$;
$\tanh(m) = \operatorname{expm1}(2m)/(\operatorname{expm1}(2m) + 2)$;
$\operatorname{asin}(m) = \operatorname{atan}(m/\sqrt{(1 - m)(1 + m)})$
and $\operatorname{acos}(m) = 2\operatorname{atan}(\sqrt{(1 - m)/(1 + m)})$
with $1 \pm m$ exact;
$\operatorname{asinh}(m) = \operatorname{log1p}(|m| + m^2/(1 + \sqrt{1 + m^2}))$;
$\operatorname{acosh}(m) = \operatorname{log1p}(d + \sqrt{d(m + 1)})$
with $d = m - 1$ exact; and
$\operatorname{atanh}(m) = \operatorname{log1p}(2|m|/(1 - |m|))/2$.
Tiny and huge arguments use the leading terms with a bounded remainder
($m$, $\log(2|m|)$, $\pm 1$), and $\operatorname{asin}(\pm 1)$,
$\operatorname{acos}(\pm 1)$ and $\operatorname{acosh}(1)$ are exact.
The radius $\rho$ enters through a bound for the variation over the
ball, for example $\rho \cosh(|m| + \rho)$ for $\sinh$ and
$\rho/(1 - |m| - \rho)$ for atanh.
:::
:::

::: function
void mp_real_exp_notab_log2(mp_real_t res, const mp_real_t x, slong n)
void mp_real_exp_notab_squaring(mp_real_t res, const mp_real_t x, slong n)
void mp_real_exp_agm(mp_real_t res, const mp_real_t x, slong n)

Set `res` to a ball containing $\exp(y)$ for every $y$ in $x$, to a
relative accuracy of about $B^{-n}$, without the precomputed tables of
`mp_real_exp_bits`. These are meant for a single exponential at huge
precision, where building the tables does not pay. The ball handling
and the range $|y| < 2^{\bits - 5}$ are those of `mp_real_exp_bits`.

::: details
- `mp_real_exp_notab_log2` reduces by $\log 2$ (computing and caching
  $\log 2$) and calls `_mp_real_exp_notab`; it is the fastest when
  $\log 2$ is already cached.
- `mp_real_exp_notab_squaring` needs no $\log 2$:
  $\exp(|m|) = \exp(|m| 2^{-s})^{2^s}$ with $|m| 2^{-s} < 1/2$, the
  reduced value from `_mp_real_exp_notab` and $s$ ball squarings (with
  $s + 16$ guard bits), and one reciprocal for $m < 0$. It is faster
  than the reduction when $\log 2$ must be computed.
- `mp_real_exp_agm` computes $y_0 = \exp(m)$ to a fraction of the
  precision by the squarings and corrects it by one Newton–Taylor step
  on the AGM logarithm of $y_0$ ([](algorithms/floating-point/elementary-functions)),
  in $O(M(n) \log n)$. It calls the squarings below about 30 limbs and
  overtakes them at about 15000 limbs.
:::
:::

## The partition function {#partitions}

::: function
void mp_real_partitions_hrr(mp_real_t res, ulong nhi, ulong nlo)
void _mp_real_partitions_hrr(mp_real_t res, const fmpz_t n, int flags)

Set `res` to the exact number of partitions $p(n)$ of
$n = n_{\mathit{hi}} B + n_{\mathit{lo}}$ (first function) or of any
`fmpz` $n$ (second function), as an exact `mp_real_t`. The size of the
result is not limited by that of an `fmpz_t`.

::: details
`mp_real_partitions_hrr` reads $p(n)$ from a table for $n < 128$, uses
the pentagonal recurrence while $p(n)$ fits a word ($n < 417$ on 64-bit
machines), and otherwise the Hardy–Ramanujan–Rademacher formula
[@Rademacher1937] [@Joh2012]

$$
p(n) = \sum_{k=1}^N A_k(n) \sqrt{\frac{3}{k}} \frac{4}{24n - 1}
\left(\cosh z_k - \frac{\sinh z_k}{z_k}\right) + R_N, \qquad
z_k = \frac{\pi\sqrt{24n - 1}}{6k},
$$

with $A_k(n)$ from `arith_hrr_expsum_factored`, $N$ and the bound for
$R_N$ from `partitions_rademacher_bound`, and each term at the precision
it contributes at. Terms of at most about 200 bits are evaluated in
batches of `dfloat` balls, the others in `mp_real` balls, with
$\exp(z_k)$ obtained as the $k$-th root of $\exp(z_1)$ (computed by
`mp_real_exp_notab_squaring` or, at huge precision, `mp_real_exp_agm`).
The resulting ball is resolved by `mp_real_unique_integer`
([](algorithms/floating-point/mp-real-balls#unique-integer)).

With several threads (`flint_set_num_threads`, from $n = 10^8$) the
constants are computed once with the whole thread budget, and the terms
are summed in ranges of $k$ of about equal cost, each into its own
accumulator. The peak memory is about that of computing $\pi$ to the
same precision.

`_mp_real_partitions_hrr` uses the formula for every $n \ge 2$; `flags`
is a combination of `MP_REAL_PARTITIONS_NO_DFLOAT` (all terms in
`mp_real`) and `MP_REAL_PARTITIONS_TEST` (the rarely taken paths at any
$n$), for testing. `partitions_fmpz_fmpz` in [](partitions) calls these
functions. The method is described in
[](algorithms/number-theory/partitions).
:::
:::

::: macro
MP_REAL_PARTITIONS_NO_DFLOAT
MP_REAL_PARTITIONS_TEST

Flags for `_mp_real_partitions_hrr`.
:::

## Hypergeometric series {#hypgeom}

The general binary splitting engine of the module, used for all the
constants. [](mp_real/hypgeom) gives examples, and
[](algorithms/floating-point/hypergeometric-series) the algorithm.

::: type
mp_real_hypgeom_int_struct

A signed integer $(-1)^{\mathit{neg}} (d, n)$ given by a limb array `d`
of length `n` (`n` = 0 for zero) and a sign `neg`.
:::

::: type
mp_real_hypgeom_series_struct

A series

$$
S = \left(\frac{c_Q + c_P \sum_{k \ge 1} \frac{P(k)}{Q(k)}
\prod_{j=1}^{k-1} \frac{R(j)}{Q(j)}}{c_D}\right)^{\pm 1},
$$

given by the fields `power` ($\pm 1$), `coefP`, `coefQ`, `coefD` (the
integers $c_P$, $c_Q$, $c_D$) and the integer polynomials `P`, `Q`, `R`
as arrays of coefficients (of degree 0 first) with lengths `Plen`,
`Qlen`, `Rlen`. This is the format of the `SeriesHypergeometric`
formula files of y-cruncher [@Yee2025].

::: details
The polynomial $Q$ must not vanish at the positive integers, $c_D \ne
0$, and the series must converge geometrically: $\deg R \le \deg Q$,
with $|\operatorname{lc}(R)| < |\operatorname{lc}(Q)|$ when the degrees
agree. Unlike y-cruncher, the coefficients may have any size. Prefactors
such as the $1/\sqrt{10005}$ of the Chudnovsky formula are left to the
caller.
:::
:::

::: function
void mp_real_hypgeom_series(mp_real_t res, const mp_real_hypgeom_series_struct * s, slong n)
void mp_real_hypgeom_series_int64(mp_real_t res, int power, int64_t coefP, int64_t coefQ, int64_t coefD, const int64_t * P, slong Plen, const int64_t * Q, slong Qlen, const int64_t * R, slong Rlen, slong n)

Set `res` to a ball containing $S$ with about $\bits(n - 1)$ accurate
bits. The second function takes the coefficients as `int64_t`, as they
appear in y-cruncher's formula files. Both throw if the series is not
geometrically convergent, if $Q(k) = 0$ for some $k \ge 1$, if $c_D =
0$, or if the tail cannot be bounded.

::: details
The series is summed by binary splitting with exact leaves of 24 to 32
terms, removal of common factors at the nodes, truncation of the values
near the root to the working precision, and transform reuse in the
merges; the number of terms comes from a rigorous bound for the tail,
which is added to the radius. Below about 30 limbs the fixed setup cost
of 1 to 2 microseconds dominates. The computation uses the threads of
`flint_set_num_threads`, and its result does not depend on their number.
:::
:::

::: function
void _mp_real_atan_frac_bsplit(mp_real_t res, nn_srcptr p, slong pn, nn_srcptr q, slong qn, slong n)
void _mp_real_atanh_frac_bsplit(mp_real_t res, nn_srcptr p, slong pn, nn_srcptr q, slong qn, slong n)

Set `res` to $\operatorname{atan}(p/q)$ resp.
$\operatorname{atanh}(p/q)$ for the integers $0 \le p < q$ given as
limb arrays, to about $n$ limbs. The quotient $p/q$ must be bounded away
from 1 (at most about 0.99); $p \ge q$ throws.

::: details
These sum the Taylor series $(p/q)\sum_k s^k (p/q)^{2k}/(2k + 1)$,
$s = \mp 1$, by binary splitting with exact leaves of up to 32 terms.
When $q^2$ is large compared to the $2k + 1$, the tree is balanced
($2^J$ leaves of $L$ terms) so that the powers $q^{2L 2^i}$ and
$p^{2L 2^i}$ needed by the merges can be tabulated by squarings, which
saves a full-size product per node. The Machin-type formulas of the
module use these functions for arguments of few bits relative to the
precision, and Zuniga's series ([](algorithms/floating-point/diophantine-reduction))
otherwise.
:::
:::

## Progress reports {#verbose}

::: function
void mp_real_set_verbose(int level)
int mp_real_get_verbose(void)

Set or read the level of progress reports on `stderr` for huge
computations: 0 is silent, 1 prints at most one line every two seconds
from the long loops and the end of every phase that took that long,
and 2 prints the start and end of every phase above a few thousand limbs
with its size, thread budget, time and (on Linux) memory use.

::: details
The default is 0 unless the environment variable
`FLINT_MP_REAL_VERBOSE` sets another level; it is read once, when no
level has been set. Each line carries the time since the first message
and the id of the calling thread; lines from different threads never
interleave.
:::
:::

## Fixed-point elementary functions {#fixed-elementary}

The kernels on fixed-point numbers that the ball functions are built
on, following [@HJ2024] and [@Joh2014c]. They take an argument in
$[0, 1)$ (or a reduced argument below $2^{-r}$) at $n$ fraction limbs
and write a bound for the error in ulps to `err`, as described under
[](#types). The guide [](mp_real/fixed-point) explains how to use them
and which kernel serves which precision; the algorithms are described in
[](algorithms/floating-point/elementary-functions).

### Series for small arguments {#series}

::: function
void _mp_real_exp_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)
void _mp_real_sin_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)
void _mp_real_cos_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)
void _mp_real_sin_cos_rs(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n)
void _mp_real_sinh_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)
void _mp_real_cosh_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)
void _mp_real_sinh_cosh_rs(nn_ptr ysinh, nn_ptr ycosh, ulong * err, nn_srcptr x, slong n)
void _mp_real_atan_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)
void _mp_real_atanh_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)

Evaluate the Taylor series of the function at $(x, n)$ with
$x < 2^{-32}$ by rectangular splitting. The outputs have $n + 1$ limbs
for exp, sin, cos, sinh and cosh and $n$ limbs for atan and atanh; the
combined versions allow either output to be `NULL`. `err` is at most 10
for exp and 15 for the others, and the bound is one-sided (the result
never exceeds the true value) for exp, sinh, cosh and atanh.

::: details
The functions dispatch on the top limbs of $x$: straight-line generated
code for $2^{-64} \le x < 2^{-32}$ and $x < 2^{-64}$ at the sizes on
the hot paths of the elementary functions (64-bit only), a generic
engine with the same calls in the same order at the other sizes, and the
tapered rectangular splitting `_mp_real_series_rs` beyond. The
precondition $x < 2^{-32}$ is checked only by `FLINT_ASSERT`.
:::
:::

### Bitwise argument reduction {#bitwise}

::: function
void _mp_real_exp_bitwise_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n, int r)
void _mp_real_log1p_bitwise_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n, int r)
void _mp_real_sin_cos_bitwise_rs(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n, int r)
void _mp_real_tan_bitwise_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n, int r)
void _mp_real_atan_bitwise_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n, int r)

Evaluate $\exp(x)$, $\log(1 + x)$, $\sin x$ and $\cos x$, $\tan x$ resp.
$\operatorname{atan} x$ for any $x \in [0, 1)$ by bitwise (BKM-type)
argument reduction below $2^{-r}$ with cached tables of
$\log(1 + 2^{-i})$ resp. $\operatorname{atan}(2^{-i})$, followed by a
series and a reconstruction. The outputs have $n + 1$ limbs for exp,
sin, cos and tan (either of `ysin`, `ycos` may be `NULL`) and $n$ limbs
for log1p and atan. `err` is at most $9r + 100$ (exp), $3r + 64$
(log1p), $6r + 128$ (sin and cos), $8r + 256$ (tan) resp. $4r + 64$
(atan), for the $r$ used. $r = 0$ selects the tuned default, at most
768; explicit values require $r \ge 32$. On 32-bit machines $n \ge 2$
is required.

::: details
All work happens at the output precision, which is why the bound grows
with $r$; a caller wanting sub-ulp accuracy pads the precision itself.
With $r = 0$, sizes $n \le 10$ on 64-bit machines go to the per-size
kernels (table-driven at one to three limbs, dense-table at four to ten
limbs), and from a few limbs the reduction runs with signed digits on
vector lanes in radix $2^{56}$ (`_mp_real_exp_bitwise_vec` and its
siblings, listed under [](#internal)). The tables are process-wide
caches built on first use at the largest precision and index range
requested so far. See [](algorithms/floating-point/elementary-functions#bitwise).
:::
:::

::: function
int _mp_real_exp_bitwise_rs_default_r(slong n)
int _mp_real_log1p_bitwise_rs_default_r(slong n)
int _mp_real_atan_bitwise_rs_default_r(slong n)
int _mp_real_trig_bitwise_rs_default_r(slong n)

Return the reduction parameter that $r = 0$ selects at $n$ limbs: the
depth of the per-size kernel where one exists, and the tuned value
beyond (from 64 at a dozen limbs up to 768). The last one serves
`_mp_real_sin_cos_bitwise_rs` and `_mp_real_tan_bitwise_rs`.
:::

### Reduced arguments {#reduced}

::: function
void _mp_real_exp_reduced(nn_ptr y, ulong * err, nn_srcptr t, slong n, flint_bitcnt_t r, int alg)
void _mp_real_sin_cos_reduced(nn_ptr ysin, nn_ptr yg, ulong * err, nn_srcptr t, slong n, flint_bitcnt_t r, int alg)

Evaluate $\exp(t)$, respectively $\sin t$ and $g = 1 - \cos t$, for a
reduced argument $(t, n)$ with $t < 2^{-r}$, $r \ge 16$ ($r \ge 32$ for
`alg` = 1, 2). The exponential is written to $n + 1$ limbs; the sine and
$g$ are pure fractions written to $n$ limbs, but both buffers need room
for $n + 1$ limbs. `err` is the rigorous bound computed at run time, at
most 96.

::: details
`alg` selects the method: 0 the tuned choice; 1 the rectangular
splitting series (`_mp_real_exp_rs`, `_mp_real_sin_cos_rs`); 2 the
$\sinh$ series and a square root, resp. the sine series and
$g = 1 - \sqrt{1 - \sin^2 t}$ (half the terms); 3 one bit-burst step,
the leading slice of $t$ by binary splitting and the rest by the series;
4 the full bit-burst algorithm, asymptotically quasi-optimal. The
cosine is returned as $g$, which keeps the information near the top for
the half-angle reconstructions that consume it. Except for the direct
series, the arithmetic around the series and the binary splitting is
`mp_real_t` ball arithmetic, so that the bound comes out rigorous. The
binary splittings share the transforms of their power tables across the
merges (see [](algorithms/floating-point/hypergeometric-series#powtab)).
:::
:::

### Table-free evaluation {#notab}

::: function
void _mp_real_exp_notab(nn_ptr y, ulong * err, nn_srcptr x, slong n)
void _mp_real_exp_notab_r(nn_ptr y, ulong * err, nn_srcptr x, slong n, int r)
void _mp_real_sin_cos_notab(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n)
void _mp_real_sin_cos_notab_r(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n, int r)

Evaluate $\exp(x)$, respectively $\sin x$ and $\cos x$, for any $(x, n)$
in $[0, 1)$ without precomputed tables: the argument is halved down to
the depth $r$, the reduced function evaluated by `_mp_real_exp_reduced`
resp. `_mp_real_sin_cos_reduced`, and the result squared resp. the angle
doubled back. Outputs have $n + 1$ limbs; `err` is at most 128. The
`_r` versions take the depth $r$ explicitly, for tuning.

::: details
With $z$ the leading zero bits of $x$, the argument is halved
$h = \max(0, r - z)$ times. The exponential is squared $h$ times (one
`sqrhigh` each, all values in $[1, e)$); relative errors double per
squaring, which $h + \log_2 n + 8$ guard bits absorb. The angle is
doubled on $g = 1 - \cos$ alone, $g(2\theta) = 2g(2 - g)$, and the sine
recovered at the end as $\sqrt{2g - g^2}$; $2h + z + \log_2 n + 8$ guard
bits absorb the amplification. The tuned depth is 16 or 32 for the
exponential and 24 or 32 for the sine and cosine.
:::
:::

### Diophantine argument reduction {#diophantine}

::: function
void _mp_real_exp_diophantine(nn_ptr y, ulong * err, nn_srcptr x, slong n)
void _mp_real_exp_diophantine_tune(nn_ptr y, ulong * err, nn_srcptr x, slong n, slong num_primes, double max_weight)
void _mp_real_sin_cos_diophantine(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n)
void _mp_real_sin_cos_diophantine_tune(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n, slong num_primes, double max_weight)
void _mp_real_tan_diophantine(nn_ptr res, ulong * err, nn_srcptr x, slong n)
void _mp_real_tan_diophantine_tune(nn_ptr res, ulong * err, nn_srcptr x, slong n, slong num_primes, double max_weight)

Evaluate $\exp(x)$, $\sin x$ and $\cos x$ (either may be `NULL`), resp.
$\tan x$ for $(x, n)$ in $[0, 1)$ by multi-prime argument reduction
[@HJ2024]: $x$ minus an integer combination of $\log p_j$ (resp. of
$\pi/2$ and the arguments of Gaussian primes), the reduced function by
`_mp_real_exp_reduced` resp. `_mp_real_sin_cos_reduced`, and a product
by a rational (resp. Gaussian rational) number. Outputs have $n + 1$
limbs; `err` is at most 3 for exp and 4 for the others. The `_tune`
versions take the number of primes (2 to 64) and the weight budget of
the reduction.

::: details
The defaults are 13 primes with a budget equal to the precision in bits
for the exponential, and 32 Gaussian primes with four times that budget
for the trigonometric functions. The logarithms and angles of the primes
are process-wide caches computed by Machin-type formulas. See
[](algorithms/floating-point/diophantine-reduction).
:::
:::

### Newton–Taylor inverses and the AGM logarithm {#newton-taylor}

::: function
void _mp_real_neglog_newton(nn_ptr y, ulong * err, nn_srcptr x, slong n)
void _mp_real_neglog_newton_tune(nn_ptr y, ulong * err, nn_srcptr x, slong n, int forward, slong N)
void _mp_real_atan_newton(nn_ptr y, ulong * err, nn_srcptr x, slong n)
void _mp_real_atan_newton_tune(nn_ptr y, ulong * err, nn_srcptr x, slong n, int forward, slong N)

Set $(y, n)$ to $-\log x$ for $(x, n)$ in $[1/2, 1)$, respectively to
$\operatorname{atan} x$ for $(x, n)$ in $[0, 1)$, by one Newton–Taylor
step over the forward function: $-\log x = t - \log(1 + w)$ with
$w = x e^t - 1$, resp. $\operatorname{atan} x = t + \operatorname{atan} w$
with $w = (x\cos t - \sin t)/(x \sin t + \cos t)$, from a starting value
$t$ correct to a fraction of the precision. `err` is computed at run
time and at most 2.

::: details
The forward functions are the diophantine exponential and sine and
cosine at $n + 1$ limbs. The `_tune` versions take the forward
algorithm (0 diophantine, 1 bitwise, 2 table-free) and the number $N$
of terms of the series in $w$ (up to 16, 0 for the default); the
starting value is computed at about $1/(N + 1)$ resp. $1/(2N)$ of the
precision, by the bitwise kernels up to 512 limbs and recursively above.
The whole costs one forward evaluation plus about a tenth of one. These
are the kernels of `mp_real_log_bits` and `mp_real_atan_bits` above 600
limbs.
:::
:::

::: function
void _mp_real_neglog_agm(nn_ptr y, ulong * err, nn_srcptr x, slong n)
void _mp_real_neglog_agm_tune(nn_ptr y, ulong * err, nn_srcptr x, slong n, slong N)

Set $(y, n)$ to $-\log x$ for $(x, n)$ in $[1/2, 1)$ by the
Sasaki–Kanada AGM formula [@SK1982]; `err` is computed at run time and
at most 2. The second function takes the number $N$ of terms of the
theta series (0 for the default, 4).

::: details
With $r = x 2^{-e}$ and $q = r^4$,
$$
-\log x = \frac{\pi/4}{\operatorname{agm}(4 r^2 S_2^2, \theta_3^2)} - e \log 2,
\qquad \theta_3 = 1 + 2\sum_{k \ge 1} q^{k^2}, \quad
S_2 = 1 + \sum_{k \ge 1} q^{k(k+1)},
$$
with $\pi/4$ and $\log 2$ from the shared caches and the shift
$e = \lceil p/(4(N + 1)^2) \rceil$ for $p$ bits. The cost is
$O(M(n) \log n)$ against the $O(M(n) \log^2 n)$ of
`_mp_real_neglog_newton`; with warm caches the two cross over near
$7 \cdot 10^5$ limbs, and against the table-free exponential the AGM is
faster from about 20000 limbs. `mp_real_exp_agm` uses it.
:::
:::

### Rational multiples of pi {#sin-cos-pi}

::: function
void _mp_real_sin_cos_pi_ui_div_ui(nn_ptr ys, nn_ptr yc, ulong * err, ulong p, ulong q, slong n)
void _mp_real_sin_cos_pi_ui_div_ui_tune(nn_ptr ys, nn_ptr yc, ulong * err, ulong p, ulong q, slong n, int alg, int r)

Set $(y_s, n + 1)$ and $(y_c, n + 1)$ (either may be `NULL`) to
$\sin(\pi p/q)$ and $\cos(\pi p/q)$ for $q \ge 1$ and $0 \le p \le q/2$.
`err` is at most 4 for the closed forms and the Chebyshev iteration and
the kernel's bound plus 2 (at most $6r + 130$, $r \le 768$) when the
sine and cosine kernel is used.

::: details
After reduction of the fraction, $q \le 2$ gives exact values, $q = 3,
4, 6$ quadratic irrationals from `_mp_real_rsqrt_ui_newton`, and $q = 5,
8, 10, 12$ closed forms with square roots from about 12 limbs (24 for
$q = 8$). For denominators up to a precision-dependent bit length (none
below 48 limbs, 12 bits at 512 limbs, 48 bits at 65536 limbs, any beyond)
the value is refined as a root of a Chebyshev polynomial by a
high-order iteration costing $O(M(n) \log q)$; otherwise the kernel of
`mp_real_sin_cos_bits` is applied to $\pi p/q$. The `_tune` version
takes the algorithm (`alg` = 0 tuned, 1 kernel, 2 Chebyshev iteration,
3 closed forms) and the order of the iteration (`r` = 0 for the default
12, at most 16). See [](algorithms/floating-point/elementary-functions#chebyshev).
:::
:::

## Fixed-point Newton inverses and roots {#newton}

The *newton* suffix marks functions whose results are not ulp-accurate:
their error bounds are relative to the output. Inputs are fixed-point
numbers $(a, a_n)$ with $a_n$ fraction limbs, and outputs written to
$n + 2$ limbs have $n$ fraction limbs and two integral limbs (the top
one possibly zero). `mpn_extras`, `arf` and the ball functions of this
module build their divisions and square roots on these functions.

::: function
void _mp_real_inv_newton_basecase(nn_ptr q, nn_srcptr a, slong an, slong n)
void _mp_real_inv_newton(nn_ptr q, nn_srcptr a, slong an, slong n)

For $a \in [1/B, 1)$ with $a_{a_n - 1} \ne 0$, set $(q, n + 2)$ to an
approximation of $1/a \in (1, B]$ with absolute error at most
$4 B^{-n}/a$.

::: details
The basecase takes the truncated reciprocal of the top
$\min(a_n, n + 1)$ limbs of $a$ by `_flint_mpn_inv_basecase`. The main
function runs a Newton iteration on middle products, each step of
second order, $T' = T + T(1 - aT)$, or third order,
$T' = T(1 - e + e^2)$ with $e = aT - 1$; the third-order step is used
for $24 \le n \le 80$ and $n \ge 1000$. See
[](algorithms/integers/division).
:::
:::

::: function
void _mp_real_div_newton_invmul(nn_ptr q, nn_srcptr b, slong bn, nn_srcptr a, slong an, slong n)
void _mp_real_div_newton(nn_ptr q, nn_srcptr b, slong bn, nn_srcptr a, slong an, slong n)

For $b \in [0, 1)$ with $b_n \ge 1$ fraction limbs and $a$ as for
`_mp_real_inv_newton`, set $(q, n + 2)$ to an approximation of $b/a$
with absolute error at most $4 B^{-n}/a$. The first function multiplies
$b$ by `_mp_real_inv_newton`; the second performs a Karp–Markstein
iteration [@KM1997].
:::

::: function
void _mp_real_rsqrt_ui_newton_basecase(nn_ptr res, ulong a, slong n)
void _mp_real_rsqrt_ui_newton(nn_ptr res, ulong a, slong n)

Set $(\mathit{res}, n)$ to the fraction limbs of an approximation of
$1/\sqrt{a}$ for a word $2 \le a < B$, with error at most $2B^{-n}$. The
main function uses the third-order steps of `_mp_real_rsqrt_newton`,
with the residual $a y^2 - 1$ exact and short.
:::

::: function
void _mp_real_rsqrt_newton_basecase(nn_ptr q, nn_srcptr a, slong an, slong n)
void _mp_real_rsqrt_newton(nn_ptr q, nn_srcptr a, slong an, slong n)

For $a \in [B^{-2}, 1)$ with one of its top two limbs nonzero, set
$(q, n + 2)$ to an approximation of $1/\sqrt{a} \in (1, B]$ with
absolute error at most $4 B^{-n}/\sqrt{a}$.

::: details
The basecase combines `flint_mpn_sqrtrem` and `flint_mpn_tdiv_qr`. The main
function performs third-order steps from $T \approx a^{-1/2}$ at about
$n/3$ limbs:
$$
u = 1 - aT^2, \qquad w = \frac{u}{2} + \frac{3u^2}{8}, \qquad q = T + Tw,
$$
with $u$ formed as a band of the product, $u^2$ from its top limbs and
$Tw$ by a middle product; the neglected term $(5/16)u^3 T$ is below the
rounding. `flint_mpn_sqrtrem` takes its reciprocal square root from
this function in its Newton regime.
:::
:::

::: function
void _mp_real_sqrt_newton_rsqrtmul(nn_ptr q, nn_srcptr a, slong an, slong n)
void _mp_real_sqrt_newton(nn_ptr q, nn_srcptr a, slong an, slong n)

For $a$ as for `_mp_real_rsqrt_newton`, set $(q, n + 2)$ to an
approximation of $\sqrt{a} \in [1/B, 1)$ (the computed value can round
up to 1) with absolute error at most $4 B^{-n}/\sqrt{a}$, which is
proportional to $1/\sqrt{a}$, not to the output. The first
function multiplies by the reciprocal square root; the second performs a
Karp–Markstein iteration.
:::

## Reduction tables {#tables}

The integer relations and Machin-type formulas behind the diophantine
argument reduction. `arb` uses the same tables
(`arb_exp_arf_log_reduction`, `arb_sin_cos_arf_atan_reduction`,
`arb_log_primes_vec_bsplit`, `arb_atan_gauss_primes_vec_bsplit`). See
[](algorithms/floating-point/diophantine-reduction).

::: macro
MP_REAL_REL_MAX
MP_REAL_LOG_PRIMES_MAX
MP_REAL_REL_TERMINATOR

The largest number of primes of a relation table (64), the same bound
for the cached logarithms of primes, and the value $-32768$ that starts
the row ending a relation table.
:::

::: type
mp_real_rel_struct

A table of integer relations
$\sum_j d_{ij}\alpha_j = \epsilon_i$ among $\alpha_j = \log p_j$ for the
first `num` primes (`gaussian` = 0), or among $\alpha_0 = \pi/2$ and
$\alpha_j = 2\arg \pi_j$ for the first nonreal Gaussian primes $\pi_j$
(`gaussian` = 1), with $|\epsilon_i|$ decreasing. The fields are `num`,
`rows`, `gaussian`, `is_static`, `primes` (`NULL` for Gaussian primes),
`weights` ($\log_2 p_j$ resp. $\log N(\pi_j)$), `d` (the `rows` $\times$
`num` coefficients followed by a row starting with
`MP_REAL_REL_TERMINATOR`), `epsilon`, `epsilon_inv` and `epsilon_min`.
:::

::: function
const mp_real_rel_struct * _mp_real_rel_table(int gaussian, slong num)
int _mp_real_rel_table_is_cached(int gaussian, slong num)

Return the relation table for the first `num` primes or Gaussian
primes, $2 \le \mathit{num} \le$ `MP_REAL_REL_MAX`, respectively whether
a call would find it ready.

::: details
Tables for `num` = 2, 4, 6, 8, 10, 12, 13, 16, 20, 24, 32, 40 and 48
are stored in the library; other sizes are generated on first use
(seconds around 48 primes) and cached per thread until `flint_cleanup`.
:::
:::

::: macro
MP_REAL_MACHIN_X_LIMBS
MP_REAL_MACHIN_MAX_NUM
MP_REAL_MACHIN_MAX_PRIMES

The number of limbs storing each argument of a Machin-type set
($128/\bits$), the largest number of terms of a set (64), and the
largest number of word-size primes used to reconstruct its coefficient
matrix by modular inversion (8).
:::

::: type
mp_real_machin_struct

A Machin-type set
$$
\log p_i = \frac{1}{\mathit{den}} \sum_j c_{ij}\operatorname{atanh}\frac{1}{x_j},
\qquad
\arg \pi_i = \frac{1}{\mathit{den}} \sum_j c_{ij}\operatorname{atan}\frac{1}{x_j},
$$
with fields `num`, `x` (the arguments, `MP_REAL_MACHIN_X_LIMBS` limbs
each), `den`, `cbits` (one more than the largest bit length of a
coefficient) and `gaussian`. The coefficient matrix is not stored but
reconstructed on first use.
:::

::: function
const mp_real_machin_struct * _mp_real_machin_table(int gaussian, slong num)
slong _mp_real_machin_table_max(int gaussian)
void _mp_real_machin_get_x(fmpz_t q, const mp_real_machin_struct * tab, slong i)
void _mp_real_machin_get_c(fmpz_t c, const mp_real_machin_struct * tab, slong i, slong j)
void _mp_real_machin_get_c_row(fmpz * row, const mp_real_machin_struct * tab, slong i)
nn_srcptr _mp_real_machin_c_row_raw(const mp_real_machin_struct * tab, slong i, slong * climbs, const unsigned char ** csign)

`_mp_real_machin_table` returns the best set for `num` values and
`_mp_real_machin_table_max` the size of the largest set. The accessors
return the argument $x_i$, the coefficient $c_{ij}$, or a row of
coefficients, the last one as the cached limbs (`climbs` limbs per
coefficient) and sign bytes.

::: details
The best set is the largest with at most `num` terms, the remaining
values being obtained by one followup series each, unless that needs
more followups than a measured limit, in which case the next larger set
is used. Sets exist for every size from 4 (3 for the Gaussian primes)
to 32, and for 40 and 48; they were found with the program at
https://github.com/fredrik-johansson/machin. The coefficients are
reconstructed by factoring each argument over the primes of the set and
inverting the exponent matrix modulo one or two primes, and cached per
thread.
:::
:::

::: macro
MP_REAL_ATAN_GAUSS_MAX

The number of Gaussian primes in `_mp_real_gaussian_primes`, 64.
:::

::: variable
const signed char _mp_real_gaussian_primes[2 * MP_REAL_ATAN_GAUSS_MAX]

The real and imaginary parts, consecutively, of the first 64 nonreal
Gaussian primes in order of norm: $1 + i$, $1 + 2i$, $2 + 3i$, ...
:::

::: function
void _mp_real_atan_gauss_vec_arb(arb_ptr res, slong num, slong prec)

Set `res` to the angles $\pi/2$ and $2\arg \pi_j$, $j < \mathit{num}$,
as `arb` balls of `prec` bits.
:::

## Multithreading {#threads}

The constants, the sets of logarithms and Gaussian-prime arguments,
`mp_real_hypgeom_series` and the bit-burst evaluations behind the
elementary functions at large precision (hence the diophantine and
table-free kernels and everything built on them) use the threads of
`flint_set_num_threads`.

- A set of independent series (the terms of a Machin-type formula, the
  logarithms of Euler's constant, the parts of $\zeta(5)$ and of the
  $\Gamma$ values) runs on up to one thread per series, costliest first.
- A binary splitting splits its thread budget between the halves of a
  node, as `flint_parallel_binary_splitting` does, while the exact
  values of the node stay within four times the working precision
  (`MP_REAL_PAR_CAP`); above, the halves run one after the other with
  the whole budget and the products of a merge run in two groups on two
  threads.
- The slices of a bit-burst evaluation run as tasks, dealt in advance to
  one lane per thread, each lane folding its factors into its own
  running product. The assignment depends only on the thread setting, so
  the result is the same on every run with the same setting.

The constants and series give identical results for any number of
threads. The elementary functions at large precision choose between a
burst step and a series by the thread count, so their rounding may
depend on it, as in `arb`.

The run-time tables of the module (the bitwise tables, the run-time
levels of the dense tables, the logarithms and angles of the
diophantine reductions, and the constants) are process-wide shared
caches (`flint_shared_cache_struct`, see [](thread_support)): one
read-only version at a time, built at the largest precision requested
so far, which every thread reads; a request beyond it builds a larger
version (at least 1.5 times as long) under the cache's lock. Since the
entries are exact floors, which nest, results do not depend on which
version served them. Superseded versions stay allocated until
`flint_cleanup_master`. With the tuned defaults, the elementary
functions up to 8 limbs and the four static constants up to 64 limbs
touch no run-time table at all.

## Internal functions {#internal}

The functions below are declared in `src/mp_real/impl.h` (or defined
inline in `src/mp_real/small.h`), which the tests, tuning and profiling
programs include as `"mp_real/impl.h"`. They are not part of the stable
interface; they are listed because the algorithm pages refer to them
and because they are useful when tuning.

### Per-size kernels {#per-size}

::: function
void _mp_real_small_exp_1_mpn(nn_ptr y, ulong * err, nn_srcptr v)
void _mp_real_small_exp_2_mpn(nn_ptr y, ulong * err, nn_srcptr v)
void _mp_real_small_exp_3_mpn(nn_ptr y, ulong * err, nn_srcptr v)
void _mp_real_small_log1p_1_mpn(nn_ptr y, ulong * err, nn_srcptr v)
void _mp_real_small_log1p_2_mpn(nn_ptr y, ulong * err, nn_srcptr v)
void _mp_real_small_log1p_3_mpn(nn_ptr y, ulong * err, nn_srcptr v)
void _mp_real_small_sin_cos_1_mpn(nn_ptr ys, nn_ptr yc, ulong * err, nn_srcptr v)
void _mp_real_small_sin_cos_2_mpn(nn_ptr ys, nn_ptr yc, ulong * err, nn_srcptr v)
void _mp_real_small_sin_cos_3_mpn(nn_ptr ys, nn_ptr yc, ulong * err, nn_srcptr v)
void _mp_real_small_atan_1_mpn(nn_ptr y, ulong * err, nn_srcptr v)
void _mp_real_small_atan_2_mpn(nn_ptr y, ulong * err, nn_srcptr v)
void _mp_real_small_atan_3_mpn(nn_ptr y, ulong * err, nn_srcptr v)

The table-driven kernels at one, two and three fraction limbs (64-bit
only): for $v \in [0, 1)$, $\exp(v)$ ($n + 1$ limbs), $\log(1 + v)$
($n$ limbs), $\sin v$ and $\cos v$ ($n + 1$ limbs each) resp.
$\operatorname{atan} v$ ($n$ limbs), with `err` between 8 and 64.

::: details
They split $v$ into 7-bit (atan: 9-bit) digits against static tables of
values at three limbs and finish with short Taylor polynomials in
registers, with tapered precision; at one limb every product is a single
high word. The scalar cores (`_mp_real_small_exp_2` and so on, on limbs
passed in registers) also serve the fast paths of `nfloat`. See
[](algorithms/floating-point/elementary-functions#small).
:::
:::

::: function
void _mp_real_small_exp_pm_2(ulong * y2, ulong * y1, ulong * y0, ulong * z2, ulong * z1, ulong * z0, ulong v1, ulong v0)
void _mp_real_small_exp_pm_3(ulong * y3, ulong * y2, ulong * y1, ulong * y0, ulong * z3, ulong * z2, ulong * z1, ulong * z0, ulong v2, ulong v1, ulong v0)

For $v \in [0, 1)$ at two resp. three fraction limbs (most significant
limb first), set $y = \exp(v)$ and $z = \exp(-v)$, each with an
integral limb, at about 1.5 times the cost of one exponential. 64-bit
only.
:::

::: function
void _mp_real_small_divq_3(ulong * w2, ulong * w1, ulong * w0, ulong n2, ulong n1, ulong n0, ulong d2, ulong d1, ulong d0)

For three-limb fractions $N < D$ with $D$ normalised (top bit set), set
$w = \lfloor N B^3/D \rfloor$ by three 3-by-2 divisions in registers.
64-bit only.
:::

::: function
void _mp_real_small_atan_ratio_2(ulong * r1, ulong * r0, ulong y1, ulong y0, ulong x1, ulong x0)
void _mp_real_small_atan_ratio_3(ulong * r2, ulong * r1, ulong * r0, ulong y2, ulong y1, ulong y0, ulong x2, ulong x1, ulong x0)

For fractions $X \ge 1/2$ and $Y \le X$ at two resp. three limbs (most
significant limb first), set the result to $\operatorname{atan}(Y/X)$
at the same number of fraction limbs, with one division instead of a
division followed by the arctangent kernel. 64-bit only.
:::

::: function
void _mp_real_exp_opt_4(nn_ptr res, nn_srcptr x)
void _mp_real_exp_opt_5(nn_ptr res, nn_srcptr x)
void _mp_real_exp_opt_6(nn_ptr res, nn_srcptr x)
void _mp_real_exp_opt_7(nn_ptr res, nn_srcptr x)
void _mp_real_exp_opt_8(nn_ptr res, nn_srcptr x)
void _mp_real_exp_opt_9(nn_ptr res, nn_srcptr x)
void _mp_real_exp_opt_10(nn_ptr res, nn_srcptr x)
void _mp_real_log1p_opt_4(nn_ptr res, nn_srcptr x)
void _mp_real_log1p_opt_5(nn_ptr res, nn_srcptr x)
void _mp_real_log1p_opt_6(nn_ptr res, nn_srcptr x)
void _mp_real_log1p_opt_7(nn_ptr res, nn_srcptr x)
void _mp_real_log1p_opt_8(nn_ptr res, nn_srcptr x)
void _mp_real_log1p_opt_9(nn_ptr res, nn_srcptr x)
void _mp_real_log1p_opt_10(nn_ptr res, nn_srcptr x)
void _mp_real_atan_opt_4(nn_ptr res, nn_srcptr x)
void _mp_real_atan_opt_5(nn_ptr res, nn_srcptr x)
void _mp_real_atan_opt_6(nn_ptr res, nn_srcptr x)
void _mp_real_atan_opt_7(nn_ptr res, nn_srcptr x)
void _mp_real_atan_opt_8(nn_ptr res, nn_srcptr x)
void _mp_real_atan_opt_9(nn_ptr res, nn_srcptr x)
void _mp_real_atan_opt_10(nn_ptr res, nn_srcptr x)
void _mp_real_sin_cos_opt_4(nn_ptr ysin, nn_ptr ycos, nn_srcptr x)
void _mp_real_sin_cos_opt_5(nn_ptr ysin, nn_ptr ycos, nn_srcptr x)
void _mp_real_sin_cos_opt_6(nn_ptr ysin, nn_ptr ycos, nn_srcptr x)
void _mp_real_sin_cos_opt_7(nn_ptr ysin, nn_ptr ycos, nn_srcptr x)
void _mp_real_sin_cos_opt_8(nn_ptr ysin, nn_ptr ycos, nn_srcptr x)
void _mp_real_sin_cos_opt_9(nn_ptr ysin, nn_ptr ycos, nn_srcptr x)
void _mp_real_sin_cos_opt_10(nn_ptr ysin, nn_ptr ycos, nn_srcptr x)

The dense-table kernels at $n = 4, \ldots, 10$ fraction limbs (64-bit
only): for $x \in [0, 1)$, $\exp(x)$ ($n + 1$ limbs),
$\log(1 + x)$ and $\operatorname{atan}(x)$ ($n$ limbs), and $\sin x$
and $\cos x$ ($n + 1$ limbs each), within the bounds of the bitwise
kernels for the reduction depths `MP_REAL_*_OPT_R` of `impl.h`.

::: details
Each level of the reduction removes seven bits of the argument with an
exact one-limb multiplier against static tables, and the remainder is
finished by a series built for its exact size. The code is generated
straight-line code (one file per size, `exp_opt_<n>.c` and so on,
emitted and tuned by `dev/tune_mp_real.py`). See
[](algorithms/floating-point/elementary-functions#dense).
:::
:::

::: function
void _mp_real_exp_dense_rt(nn_ptr y, nn_srcptr v, slong n)
void _mp_real_log1p_dense_rt(nn_ptr y, nn_srcptr v, slong n)
void _mp_real_atan_dense_rt(nn_ptr y, nn_srcptr v, slong n)
void _mp_real_sin_cos_dense_rt(nn_ptr ysin, nn_ptr ycos, nn_srcptr v, slong n)

The dense-table kernels as generic loops, for $2 \le n \le$
`MP_REAL_DENSE_RT_MAX` = 32 fraction limbs (64-bit only), with the
deeper table levels computed at run time as they are needed and shared
by all threads. Outputs as for the per-size kernels. The dispatchers
call them from eleven limbs on.
:::

::: function
slong _mp_real_dense_rt_prepare(int which, slong n)

Fill every entry not yet computed of the run-time dense tables of a
family (`MP_REAL_DENSE_RT_EXP`, `MP_REAL_DENSE_RT_LOG1P` or
`MP_REAL_DENSE_RT_ATAN`, the last shared with sine and cosine) in the
tier serving $n$ limbs, and return the width of that tier: 21 for
$n \le 21$, else `MP_REAL_DENSE_RT_MAX`. 64-bit only.
:::

::: function
void _mp_real_tan_sc(nn_ptr y, ulong * err, nn_srcptr v, slong n)

Set $(y, n + 1)$ to $\tan v$ for $v \in [0, 1)$ and $n \le 32$ from
the sine and cosine kernels and one division (64-bit only). With both
kernel outputs within $E$ ulps, `err` is $\lfloor 19E/4 \rfloor + 2$.
The tangent kernel uses it up to 22 limbs.
:::

### Series {#internal-series}

::: macro
MP_REAL_SERIES_TAN
MP_REAL_SERIES_ATAN
MP_REAL_SERIES_ATANH
MP_REAL_SERIES_SIN
MP_REAL_SERIES_COS
MP_REAL_SERIES_SINH
MP_REAL_SERIES_COSH
MP_REAL_SERIES_EXP
MP_REAL_SERIES_EXP_NEG

The function selectors of `_mp_real_series_rs` and
`_mp_real_series_tapered`. `COS` and `COSH` select $1 - \cos x$ and
$\cosh x - 1$.
:::

::: function
void _mp_real_series_rs(nn_ptr res, nn_srcptr x, slong n, flint_bitcnt_t r, int func)
void _mp_real_series_rs_sin_cos(nn_ptr ysin, nn_ptr yg, nn_srcptr x, slong n, flint_bitcnt_t r, int hyperbolic)

Tapered rectangular splitting for a fraction $0 \le x < 2^{-r}$,
$r \ge 8$, of $n$ limbs. The first function computes $\exp(\pm x)$
($n + 1$ limbs, within 5 ulps) or $\sin x$, $\sinh x$,
$\operatorname{atan} x$, $\operatorname{atanh} x$, $1 - \cos x$,
$\cosh x - 1$ ($n$ limbs, within 2 ulps) as selected by `func`. The
second computes the sine and $1 - \cos$ (or $\sinh$ and $\cosh - 1$
when `hyperbolic` is set) sharing the powers; either output may be
`NULL`. Outputs may alias $x$.

::: details
With $v = x$ (exp) or $v = x^2$, the series $\sum_{k<N} a_k v^k$ is
split into blocks of $m \approx \sqrt{N/1.5}$ terms, the powers
$v, \ldots, v^m$ computed once and every term a single `mpn_addmul_1`
or `mpn_submul_1` by a one-limb coefficient. Block $b$ enters the sum
multiplied by $v^{mb} < 2^{-smb}$ ($s = r$ resp. $2r$), so it runs in a
window shortened by $smb$ bits. See
[](algorithms/floating-point/elementary-functions#series).
:::
:::

::: function
void _mp_real_series_tapered(nn_ptr res, nn_srcptr x, slong n, flint_bitcnt_t r, int func)

Set $(\mathit{res}, n)$ to $\tan x$, $\operatorname{atan} x$,
$\operatorname{atanh} x$, $\sin x$ or $1 - \cos x$ (selected by `func`)
for $0 \le x < 2^{-r}$, within 6 ulps, by Horner's rule in $x^2$ over
tapered widths with static coefficient tables. The supported ranges are
$n \le 80$, $r \ge 32$ for the tangent, $n \le 13$, $r \ge 10$ for atan
and atanh, and $n \le 16$, $r \ge 18$ for the sine and cosine.
:::

::: function
void _mp_real_series_rs_tan(nn_ptr res, nn_srcptr x, slong n, flint_bitcnt_t r)
int _mp_real_series_rs_tan_ok(slong n, flint_bitcnt_t r)

Set $(\mathit{res}, n)$ to $\tan x$ for $0 \le x < 2^{-r}$, $r \ge 32$,
within 3 ulps, for $n \le 80$ and beyond as long as
`_mp_real_series_rs_tan_ok` returns nonzero. `res` may alias $x$.

::: details
The coefficients $c_k = T_{k+2}/(2k+3)!$ of $\tan t = t + t^3\sum_k c_k
t^{2k}$ ($T_m$ the tangent numbers) are not ratios of small integers,
but the leading ones share small denominators: chunk $j$ holds the terms
whose common denominator $Q_j$ fits $j + 1$ limbs (12 terms in one limb,
then 5 to 7 more per limb on 64-bit machines), with integer numerators
in static tables. One tapered rectangular splitting runs over the
chunks, carrying $Q_j$ times the partial sum; crossing into chunk
$j - 1$ costs one division by $Q_j/Q_{j-1}$. The tail is the tapered
Horner scheme of `_mp_real_series_tapered`.
:::
:::

### Vector reductions {#vector}

::: function
void _mp_real_exp_bitwise_vec(nn_ptr res, nn_srcptr x, slong n, int r)
void _mp_real_log1p_bitwise_vec(nn_ptr res, nn_srcptr x, slong n, int r)
void _mp_real_atan_bitwise_vec(nn_ptr res, nn_srcptr x, slong n, int r)
void _mp_real_tan_halfangle_vec(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan, nn_srcptr x, slong n, int r)

The bitwise reductions of `_mp_real_exp_bitwise_rs`,
`_mp_real_log1p_bitwise_rs`, `_mp_real_atan_bitwise_rs` and of the
sine, cosine and tangent, with signed digits and the full-precision
arithmetic in radix $2^{56}$ on `vec4n` lanes, for $x \in [0, 1)$ and
$32 \le r \le 64n - 16$, within the bounds of the corresponding
`_bitwise_rs` functions (64-bit only). Outputs: $(\mathit{res}, n + 1)$
for exp, $(\mathit{res}, n)$ for log1p and atan, and $n + 1$ limbs for
each of $\sin$, $\cos$, $\tan$ (each may be `NULL`).

::: details
The `_bitwise_rs` functions call these from
`MP_REAL_EXP_BITWISE_VEC_MIN` = 4, `MP_REAL_LOG1P_BITWISE_VEC_MIN` =
48, `MP_REAL_ATAN_BITWISE_VEC_MIN` = 4 and
`MP_REAL_TRIG_BITWISE_VEC_MIN` = 4 limbs (beyond the per-size kernels
for $r = 0$). See
[](algorithms/floating-point/elementary-functions#signed-digits).
:::
:::

### Tables and series of the reductions {#internal-tables}

::: function
void _mp_real_log1p_2mexp_ui_bs(nn_ptr res, ulong i, slong n)
void _mp_real_atan_2mexp_ui_bs(nn_ptr res, ulong i, slong n)

Set $(\mathit{res}, n)$ to the floor of $\log(1 + 2^{-i})$ resp.
$\operatorname{atan}(2^{-i})$ scaled by $B^n$, or one less, for
$i \ge 1$, by binary splitting in ball arithmetic. These compute
entries of the bitwise tables.
:::

::: function
void _mp_real_log_primes_vec(mp_real_struct * res, slong num, slong n)
void _mp_real_atan_gauss_vec(mp_real_struct * res, slong num, slong n)

Set `res` to the balls $\log p_j$, resp. $\pi/2$ followed by
$2\arg \pi_j$, $j < \mathit{num}$, to about $n$ limbs, by the
Machin-type sets of `_mp_real_machin_table` and followup series, the
series running on the available threads.
:::

::: function
void _mp_real_log_ratio_zuniga(mp_real_t res, const fmpz_t u, const fmpz_t v, slong n)
void _mp_real_atan_ratio_zuniga(mp_real_t res, const fmpz_t d, const fmpz_t t, slong n)

Set `res` to $\log(u/v)$ for integers $u > v > 0$, resp. to
$\operatorname{atan}(d/t)$ for integers $0 < d < t$, by Zuniga's Ramanujan-type series
[@Zun2025] [@Zun2025b], to about $n$ limbs. See
[](algorithms/floating-point/diophantine-reduction#zuniga).
:::

::: function
void _mp_real_exp_sum_bs_powtab(nn_ptr T, slong * tn, nn_ptr Q, slong * qn, slong * QE, nn_srcptr xp, slong xn, slong D, slong N)
void _mp_real_sin_cos_sum_bs_powtab(nn_ptr A, slong * an, slong * ae, nn_ptr B, slong * bn, slong * be, nn_ptr Q, slong * qn, slong * QE, nn_srcptr xp, slong xn, slong D, slong N, slong lmax)
slong _mp_real_exp_bs_num_terms(flint_bitcnt_t r, slong prec)

The exact binary splittings of the bit-burst steps. The first computes
$T/(Q B^{QE}) = \sum_{k=1}^N x^k/(k!\,B^{Dk})$ for the integer $x$ given
by $(x_p, x_n)$, $x_n \le D$; the second the sine and $1 - \cos$ of
$x B^{-D}$ over one common denominator,
$1 - \cos = A B^{ae}/(Q B^{QE})$ and
$\sin = x B^{-D}(Q B^{QE} - B B^{be})/(Q B^{QE})$, with $A$ and $B$
truncated to `lmax` limbs. The third returns the number of terms for a
reduced argument below $2^{-r}$ and `prec` bits.

::: details
Both splittings use a table of the powers $x^{s}$ (resp. $(x^2)^s$)
for the split steps $s$ of the fixed-midpoint recursion, and compute the
`fft_small` transform of each table entry used by at least two merges
once per tree (`_mp_real_powtab_pre_init`). The output buffers must
have room as stated in `exp_sum_bs.c` and `sin_cos_sum_bs.c`.
:::
:::

::: type
mp_real_powtab_pre_struct

Transformed entries of a power table: the fields `pre` (an array of
`fft_small_mpn_mul_precomp_struct`, or `NULL` when no entry is
transformed), `storage` (one block holding all the transforms) and
`len` (the number of entries).
:::

::: function
void _mp_real_bs_step_uses(slong * uses, slong * rmax, const slong * xexp, slong len, slong N)

For a binary splitting of $N$ terms with fixed midpoints (ranges of one
or two terms being leaves) whose split steps are the `len` entries of
`xexp`, set `uses[i]` to the number of merges using the step `xexp[i]`
and `rmax[i]` to the longest right child of those merges.

::: details
The ranges of a level differ in length by at most one, so the counts
follow from $N$ alone by tracking at most two distinct lengths per
level, in $O(\log N)$ time.
:::
:::

::: function
void _mp_real_powtab_pre_init(mp_real_powtab_pre_struct * T, nn_srcptr * xpow, const slong * xlen, const slong * uses, const slong * ymax, slong len)
void _mp_real_powtab_pre_clear(mp_real_powtab_pre_struct * T)
int _mp_real_powtab_pre_mul(nn_ptr z, const mp_real_powtab_pre_struct * T, slong i, nn_srcptr y, slong yn)

Prepare, free, and use the transforms of the power-table entries
`xpow[i]` of `xlen[i]` limbs for products by partners of up to
`ymax[i]` limbs. `_mp_real_powtab_pre_mul` sets $z$ to
$\mathit{xpow}[i] \cdot y$ ($\mathit{xlen}[i] + y_n$ limbs, $z$ may
overlap $y$) through the transform of entry $i$ and returns 1, or
returns 0 without writing anything when the entry has no transform or
the product is outside its range.

::: details
An entry is transformed only when at least two merges use it
(`uses[i]` $\ge 2$), its products are in the FFT range (with the
threshold `MP_REAL_POWTAB_PRE_CUTOFF` = 150 limbs), and they have at
most `MP_REAL_POWTAB_PRE_MAX_LIMBS` = $2^{20}$ limbs, which bounds the
extra memory to some 50 MB per tree. Without `fft_small` the functions
do nothing and `_mp_real_powtab_pre_mul` always returns 0. The
transforms are read only, so the two halves of a merge may use them on
two threads. See
[](algorithms/floating-point/hypergeometric-series#powtab).
:::
:::

::: macro
MP_REAL_POWTAB_PRE_CUTOFF
MP_REAL_POWTAB_PRE_MAX_LIMBS

The FFT-range threshold (150 limbs) and the size limit ($2^{20}$ limbs
per product) of the transformed power tables. Both can be overridden at
compile time.
:::
