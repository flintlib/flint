---
title: radix.h
header: radix.h
description: Multiprecision arithmetic in a general radix such as $10^{19}$ or a power of a prime, with conversion to and from binary.
category: Integers
order: 7
---

This module implements multiprecision natural number and integer
arithmetic in a general radix $B < \beta$, where $\beta = 2^{\bits}$ is
the machine word radix. The radix is a power $B = b^e$ of a *digit radix*
$b \ge 2$, so that one limb holds $e$ digits; for decimal arithmetic one
normally takes $B = 10^{19}$ on 64-bit machines ($10^9$ on 32-bit
machines), the largest power of 10 below $\beta$.

Arithmetic in radix $B$ is somewhat slower than in the radix $\beta$, but
everything related to the base $b$ becomes cheap: reduction modulo $b^n$,
truncation to $n$ digits, shifts by digits and conversion to a string in
base $b$ take time $O(n)$ instead of $O(M(n))$ or $O(M(n) \log n)$. The two
principal applications are decimal arithmetic ($b = 10$, see the
[`decimal`](decimal) module) and $p$-adic arithmetic ($b = p$ a small
prime, with Hensel division and $p$-adic square roots; see the
[`padic_radix`](padic_radix) module).

The module provides low-level functions on limb arrays (in the style of
`mpn_extras`), conversion between radix $B$ and the binary radix
$\beta$ in quasi-linear time, which `flint_mpn_get_str` uses for decimal
output, and a memory-managed signed integer type `radix_integer_t` with a
generic ring (`gr`) context. The functions are not optimised for $b$ a
power of two, where binary arithmetic should be used. The guide
[](radix/getting-started) shows how to use the module; [](radix/conversion) describes the conversions to and from binary.

## Radix objects {#radix}

::: type
radix_struct
radix_t

A radix context, holding the digit radix $b$ and the limb radix $B = b^e$
as `nmod_t` values (fields `b` and `B`), the exponent `exp` $= e$, and
precomputed data: the powers $b^k$ for $0 \le k \le e$, data for division
by these powers (`n_div_precomp_t`), data for division by $e$ (see
`radix_div_exp`), the 2-adic valuation of $b$, inverses of the powers of
the odd part of $b$, and a function computing digit valuations. A `radix_t` is passed by reference to every function of the
module.
:::

::: macro
DIGIT_RADIX(radix)
LIMB_RADIX(radix)

The digit radix $b$ and the limb radix $B$ of a `radix_t`, as `ulong`.
:::

::: function
ulong radix_div_exp(ulong x, const radix_t radix)
ulong radix_divrem_exp(ulong * r, ulong x, const radix_t radix)

Return $\lfloor x / e \rfloor$ (and set $r = x \bmod e$), converting a
number of digits to a number of limbs without a hardware division. For
$e > 2$ this is one multiplication by a precomputed
$\lceil 2^{\bits - 2 + \ell} / e \rceil$, $\ell = \lceil \log_2 e \rceil$, and a
shift, exact for $x < 2^{\bits - 2}$; larger $x$ use `n_div_precomp`.
:::

::: function
void radix_init(radix_t radix, ulong b, unsigned int exp)

Initialises `radix` with digit radix $b$ and limb radix $B = b^e$ where
$e$ = `exp`. Requires $b \ge 2$, $e < \bits$ and $b^e < \beta$, and throws
an exception otherwise. With $e = 0$, the largest $e$ with $b^e < \beta$
is chosen.
:::

::: function
void radix_clear(radix_t radix)

Frees the memory used by `radix`.
:::

::: function
void radix_init_randtest(radix_t radix, flint_rand_t state)
void radix_init_randtest_prime(radix_t radix, flint_rand_t state)

Initialise `radix` with a random digit radix $b$ (a random prime $b$ for
the second function) and a random exponent, for testing.
:::

::: function
ulong radix_digit_radix(const radix_t radix)
ulong radix_limb_radix(const radix_t radix)
ulong radix_limb_exponent(const radix_t radix)

Return $b$, $B$ and $e$.
:::

## Natural numbers {#natural}

The low-level functions work on arrays of limbs (`ulong`), with the same
notation as in [`mpn_extras`](mpn_extras): $(x, n)$ is the natural number
$\sum_{i<n} x_i B^i$ in $[0, B^n - 1]$. Unless stated otherwise:

- input limbs must be reduced to $[0, B - 1]$, and output limbs are
  reduced;
- aliasing between outputs and inputs is allowed only if the aliased
  operands start at the same address;
- *normalised* means that the top limb is nonzero.

::: function
ulong _radix_size_digits_1(ulong c, const radix_t radix)
ulong radix_size_digits_1(ulong c, const radix_t radix)

Return the number of digits of the limb $c$ (the least $k$ with
$c < b^k$). The underscore version requires $c > 0$; the other returns 0
for $c = 0$.
:::

::: function
ulong _radix_size_digits(nn_srcptr x, slong n, const radix_t radix)
ulong radix_size_digits(nn_srcptr x, slong n, const radix_t radix)

Return the number of digits of $(x, n)$. The underscore version requires
$n > 0$ and $x_{n-1} \ne 0$; the other allows $n = 0$ and zero top limbs,
which it skips, and returns 0 for zero.
:::

::: function
ulong _radix_valuation_digits_1(ulong c, const radix_t radix)

Returns the number of trailing zero digits of the limb $c > 0$, that is,
the largest $k$ such that $b^k$ divides $c$.
:::

::: function
void radix_rand_limbs(nn_ptr res, flint_rand_t state, slong n, const radix_t radix)
void radix_rand_digits(nn_ptr res, flint_rand_t state, slong n, const radix_t radix)

Set $(\mathit{res}, n)$ to a uniformly random integer in $[0, B^n - 1]$,
respectively $(\mathit{res}, \lceil n/e \rceil)$ to a uniformly random
integer in $[0, b^n - 1]$.
:::

::: function
void radix_randtest_limbs(nn_ptr res, flint_rand_t state, slong n, const radix_t radix)
void radix_randtest_digits(nn_ptr res, flint_rand_t state, slong n, const radix_t radix)

Like `radix_rand_limbs` and `radix_rand_digits`, but with a
non-uniform distribution producing long runs of limbs equal to 0 or
$B - 1$, for testing.
:::

## Addition and shifts {#addition}

::: function
ulong radix_add(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix)
ulong radix_sub(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix)

Set $(\mathit{res}, \mathit{an})$ to $a + b$, respectively $a - b$, modulo
$B^{an}$ and return the carry, respectively the borrow (0 or 1). Requires
$\mathit{an} \ge \mathit{bn} \ge 1$.
:::

::: function
ulong radix_neg(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)

Sets $(\mathit{res}, \mathit{an})$ to the $B$-complement $(-a) \bmod
B^{an}$ and returns the borrow (1 if $a \ne 0$, 0 otherwise).
:::

::: function
ulong radix_lshift_digits(nn_ptr res, nn_srcptr a, slong n, unsigned int e, const radix_t radix)
ulong radix_rshift_digits(nn_ptr res, nn_srcptr a, slong n, unsigned int e, const radix_t radix)

Shift $(a, n)$ by $k$ = `e` digits, for $1 \le k$ less than the limb
exponent: the left shift sets $(\mathit{res}, n)$ to $a b^k \bmod B^n$ and
returns the $k$ digits shifted out (as the low digits of the returned
limb); the right shift sets $(\mathit{res}, n)$ to $\lfloor a / b^k
\rfloor$ and returns the $k$ digits shifted out as the *high* digits of the
returned limb, that is, $(a \bmod b^k) b^{e - k}$. These are faster than
`radix_mul_1` and `radix_divrem_1` with $b^k$. Shifts by whole limbs are
pointer offsets.
:::

::: function
ulong radix_addlsh(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, unsigned int sh, const radix_t radix)
ulong radix_sublsh(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, unsigned int sh, const radix_t radix)
ulong radix_lshsub(nn_ptr res, nn_srcptr b, slong bn, unsigned int sh, nn_srcptr a, slong an, const radix_t radix)

Fused addition and subtraction with a shift of $b$ by $1 \le
\mathit{sh} < e$ digits: set $(\mathit{res}, L)$ with $L = \max(\mathit{an},
\mathit{bn} + 1)$ to $a + b \cdot b^{\mathit{sh}}$, $a - b \cdot
b^{\mathit{sh}}$ and $b \cdot b^{\mathit{sh}} - a$, and return the carry
(0 or 1), respectively the borrow, which is 0 since the result must be
nonnegative. Requires $\mathit{bn} \ge 1$. `res` may alias $a$ but not $b$.
The shifted operand is never formed in memory.
:::

::: function
int radix_cmplsh(nn_srcptr a, slong an, nn_srcptr b, slong bn, unsigned int sh, const radix_t radix)

Returns the sign of $a - b \cdot b^{\mathit{sh}}$ for normalised
$(a, \mathit{an})$ and $(b, \mathit{bn})$ and $1 \le \mathit{sh} < e$.
:::

## Multiplication {#multiplication}

::: function
ulong radix_mul_1(nn_ptr res, nn_srcptr a, slong n, ulong c, const radix_t radix)

Sets $(\mathit{res}, n)$ to $ac \bmod B^n$ and returns the carry
$\lfloor ac / B^n \rfloor$. Requires $n \ge 1$ and $c < B$.
:::

::: function
ulong radix_mul_two(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)

Sets $(\mathit{res}, \mathit{an})$ to $2a \bmod B^{an}$ and returns the
carry.
:::

::: function
void radix_mul(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix)
void radix_sqr(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)

Set $(\mathit{res}, \mathit{an} + \mathit{bn})$ to $ab$, respectively
$(\mathit{res}, 2\mathit{an})$ to $a^2$. Requires $\mathit{an},
\mathit{bn} \ge 1$; `res` must not overlap the inputs. Both are calls to
`radix_mulmid` with the full window (there is no special squaring code
yet).
:::

::: function
void radix_mulmid(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix)

Windowed product: for $\mathit{an}, \mathit{bn} \ge 1$ and
$0 \le \mathit{lo} < \mathit{hi} \le \mathit{an} + \mathit{bn}$, sets
$(\mathit{res}, \mathit{hi} - \mathit{lo})$ to the limbs
$[\mathit{lo}, \mathit{hi})$ of the product, computed without the carries
from the limbs below $\mathit{lo}$. `res` must not overlap the inputs.

::: details
Viewing $a$ and $b$ as polynomials $A, Y \in \Z[T]$ evaluated at $T = B$,
the function computes the slice
$Z_{lo:hi} = \sum_{k=lo}^{hi-1} ([T^k] AY)\, T^{k-lo}$ of the polynomial
product and writes $Z_{lo:hi}(B) \bmod B^{hi-lo}$: carries from lower
terms are omitted and carries beyond the top are discarded. In particular:

- for $\mathit{lo} = 0$ this is the low product $ab \bmod B^{hi}$;
- for $\mathit{hi} = \mathit{an} + \mathit{bn}$ this is an approximate
  high product $r$ with
  $ab \ge r B^{lo} \ge ab - \min(\mathit{an}, \mathit{bn}, \mathit{lo})
  B^{lo+1}$;
- for $\mathit{lo} = 0$ and $\mathit{hi} = \mathit{an} + \mathit{bn}$ it
  is the full product.

The classical algorithm is used if either operand or the window is
shorter than 80 limbs, and otherwise if the number of limb products the
window needs (counted by `_radix_mulmid_pairs_below`) is less than
`RADIX_MULMID_FFT_CUTOFF_RATIO` $\cdot L$; the FFT (or Kronecker
substitution without `fft_small`) is used otherwise. The cost of the FFT
depends essentially on the length $L$ of the convolution it computes,
while the classical algorithm benefits from short windows: a high
product needs about half as many limb products as the full product.
With `fft_small`, $L$ is `_radix_mulmid_fft_length`$(\mathit{an},
\mathit{bn}, \mathit{lo}, \mathit{hi})$, which is smaller than
$\mathit{an} + \mathit{bn}$ for windows that do not need all limbs or
that are computed as transposed products; without it,
$L = \mathit{an} + \mathit{bn}$.
:::
:::

::: function
void radix_mulmid_classical(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix)
void radix_mulmid_fft_small(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix)
void radix_mulmid_KS(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix)
void radix_mulmid_naive(nn_ptr res, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong lo, slong hi, const radix_t radix)

The algorithms behind `radix_mulmid`, with its interface: the schoolbook
algorithm with delayed reduction of the column sums, efficient when an
operand or the window is short; multiplication of the limb sequences as
polynomials with the small-prime FFT of `fft_small` (only available if
FLINT is built with it), reconstructing only the coefficients that reach
the window and propagating carries in radix $B$; Kronecker substitution
into one product in the binary radix $\beta$ (each limb padded to two
words); and an unoptimised reference implementation through `fmpz_poly`.

::: details
`radix_mulmid_fft_small` first drops the limbs of each operand above
$\mathit{hi}$ and the low limbs that meet the other operand only below
$\mathit{lo}$ (they contribute only carries from below the window,
which are omitted anyway). If the remaining window is shorter than the
longer operand, it is computed as a transposed product of the reversed
shorter operand and the longer one, at the cost of a product of the
shorter operand and the window, provided this transform is strictly
shorter than that of the window of the full product (not for squarings).
The middle products of the Newton iterations of the division and
inversion code (an operand of about $2n$ limbs, a window of about $n$
limbs) are the typical case. The carries are propagated from
$\mathit{lo}$ in both cases, so the result is the same. See
[](algorithms/polynomials/middle-product#transposed-fft).
:::
:::

::: function
double _radix_mulmid_pairs_below(slong an, slong bn, slong k)

Returns the number of pairs $(i, j)$ with $0 \le i < \mathit{an}$,
$0 \le j < \mathit{bn}$ and $i + j < k$, as a `double` (it may exceed the
range of `slong` on 32-bit machines). It counts the limb products of the
classical algorithm and is used for tuning decisions.
:::

::: macro
RADIX_MULMID_FFT_CUTOFF_RATIO

The tuning ratio of `radix_mulmid`. With `fft_small` it is
`FLINT_RADIX_MULMID_FFT_SMALL_CUTOFF_RATIO`, which every `flint-mparam.h`
must define: 50 on Zen 3 and Zen 4, 60 on Ice Lake and later Intel
cores and in the generic files, 80 on Skylake and Broadwell, 85 on arm64
(measured on Apple M4). The measured crossover is the same for balanced
full products, Newton middle products and short windows of a long
operand; a balanced full product of $n$ limbs switches to the FFT at
$n = 2 \cdot$ ratio, 100 to 170 limbs. Without `fft_small` the ratio is
175 (Kronecker substitution pads each limb to two words, so the
crossover for a balanced full product is at about 350 limbs).
:::

::: function
slong _radix_mulmid_fft_length(slong an, slong bn, slong lo, slong hi)

The length of the convolution that `radix_mulmid_fft_small` computes for
the window $[\mathit{lo}, \mathit{hi})$, used by `radix_mulmid` to
choose the algorithm: with $\mathit{an}'$ and $\mathit{bn}'$ the
numbers of limbs of the operands that reach the window, it is
$\min(\mathit{an}', \mathit{bn}') + (\mathit{hi}' - \mathit{lo}')$ if
the window (clipped to the product) is shorter than
$\max(\mathit{an}', \mathit{bn}')$, and $\mathit{an}' + \mathit{bn}'$
otherwise; 0 for an empty window.
:::

::: function
void _radix_mulhigh_known_low(nn_ptr out, nn_srcptr x, slong xn, nn_srcptr y, slong yn, nn_srcptr kl, slong kl_len, slong klo, slong khi, nn_ptr scratch, const radix_t radix)

The radix $B$ analogue of `_flint_mpn_mulhigh_known_low`: sets
$(\mathit{out}, \mathit{khi} - \mathit{klo})$ to the exact limbs
$[\mathit{klo}, \mathit{khi})$ of $xy$, given that its low $\mathit{klo}$
limbs are $(\mathit{kl}, \mathit{kl\_len})$ (zero above `kl_len`), using
a windowed product with guard limbs. `scratch` must have room for
$\mathit{khi}$ limbs. Used by the Hensel division and $p$-adic square
roots.
:::

## Division {#division}

Division in radix $B$ follows the structure of the binary case: small
divisors are handled directly, large ones by Newton iteration on
fixed-point approximations of the reciprocal, and medium sizes by
converting to binary and using `flint_mpn_tdiv_qr`.

::: function
ulong radix_divrem_1(nn_ptr res, nn_srcptr a, slong an, ulong d, const radix_t radix)
void radix_divexact_1(nn_ptr res, nn_srcptr a, slong an, ulong d, const radix_t radix)
ulong radix_divrem_two(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)

Set $(\mathit{res}, \mathit{an})$ to $\lfloor a/d \rfloor$ and return the
remainder, for $1 \le d < B$; the same assuming that $d$ divides $a$
(the result is undefined otherwise); and division by 2.
:::

::: function
void radix_inv_approx_basecase(nn_ptr q, nn_srcptr a, slong an, slong n, const radix_t radix)
void radix_inv_approx(nn_ptr q, nn_srcptr a, slong an, slong n, const radix_t radix)

Given $(a, \mathit{an})$ with $a_{an-1} \ne 0$, read as a fixed-point
number $a \in [1/B, 1)$ with $\mathit{an}$ fraction limbs, set
$(q, n + 2)$ to an approximation of $1/a \in (1, B]$ with $n$ fraction
limbs and two integral limbs (the top one may be zero). The absolute error
is less than $4 B^{-n}/a$ (relative error less than $4 B^{-n}$). The
general version uses Newton iteration from the basecase.
:::

::: function
void radix_div_approx_invmul(nn_ptr Q, nn_srcptr B, slong Bn, nn_srcptr A, slong An, slong n, const radix_t radix)
void radix_div_approx(nn_ptr Q, nn_srcptr B, slong Bn, nn_srcptr A, slong An, slong n, const radix_t radix)

Given the divisor $(A, \mathit{An})$ with nonzero top limb, read as a
fixed-point number in $[1/B, 1)$ with $\mathit{An} \ge 1$ fraction limbs,
and the dividend $(B, \mathit{Bn})$, read as a fixed-point number in
$[0, 1)$ with $\mathit{Bn} \ge 1$ fraction limbs, set $(Q, n + 2)$ to an
approximation of the quotient with $n$ fraction limbs and two integral
limbs; the absolute error is less than $4 B^{-n}$ divided by the divisor.
Note the order of the arguments: the dividend comes first. The
`invmul` version multiplies by an approximate reciprocal; the general
version uses Karp–Markstein iteration [@KM1997].
:::

::: function
void radix_divrem(nn_ptr q, nn_ptr r, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix)

Sets $(q, \mathit{an} - \mathit{bn} + 1)$ to $\lfloor a/b \rfloor$ and,
if $r$ is not `NULL`, $(r, \mathit{bn})$ to the remainder. Requires
$\mathit{an} \ge \mathit{bn} \ge 1$ and $b_{bn-1} \ne 0$.

::: details
One-limb divisors use `radix_divrem_1`. Divisors of at least 70 limbs,
and divisors of at least 20 limbs with quotients of at most 150 limbs,
use `radix_divrem_newton_karp_markstein`. Otherwise long quotients
($\mathit{an} > 2\mathit{bn}$) are computed blockwise by
$2\mathit{bn} \times \mathit{bn}$ divisions, and the remaining divisions
convert to binary and use `flint_mpn_tdiv_qr` (`radix_divrem_via_mpn`).
:::
:::

::: function
void radix_divrem_via_mpn(nn_ptr q, nn_ptr r, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix)
void radix_divrem_newton(nn_ptr q, nn_ptr r, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix)
void radix_divrem_newton_karp_markstein(nn_ptr Q, nn_ptr R, nn_srcptr A, slong An, nn_srcptr B, slong Bn, const radix_t radix)

The algorithms behind `radix_divrem`, with its interface: conversion to
binary, division by a Newton reciprocal, and Karp–Markstein division
(quotient from a half-precision reciprocal and one correction step).
:::

::: function
void radix_divrem_preinv(nn_ptr q, nn_ptr r, nn_srcptr a, slong an, nn_srcptr b, slong bn, nn_srcptr binv, slong binvn, const radix_t radix)

`radix_divrem` with a precomputed inverse of $b$ given as
$(\mathit{binv}, \mathit{binvn} + 2)$ with $\mathit{binvn}$ fraction limbs,
as computed by `radix_inv_approx`. Requires $\mathit{binvn} \ge
\mathit{an} - \mathit{bn} + 1$; a few extra limbs can improve performance.
$r$ may be `NULL`.
:::

::: function
int radix_div(nn_ptr q, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix)

If $b$ divides $a$, sets $(q, \mathit{an} - \mathit{bn} + 1)$ to the
quotient and returns 1; otherwise returns 0, the contents of $q$ being
undefined.
:::

::: function
void radix_divexact(nn_ptr q, nn_srcptr a, slong an, nn_srcptr b, slong bn, const radix_t radix)

Sets $(q, \mathit{an} - \mathit{bn} + 1)$ to $a/b$, assuming that the
division is exact; the result is undefined otherwise.
:::

::: function
int radix_cmp_bn_half(nn_srcptr x, slong n, const radix_t radix)

Returns $-1$, 0 or 1 according to whether $(x, n)$ is less than, equal
to or greater than $\lfloor B^n/2 \rfloor$. Requires $n \ge 1$. Used for
rounding to nearest and for symmetric remainders.
:::

## Square roots {#sqrt}

::: function
void radix_sqrtrem(nn_ptr s, nn_ptr r, nn_srcptr a, slong an, const radix_t radix)

Sets $(s, \lceil \mathit{an}/2 \rceil)$ to $\lfloor \sqrt a \rfloor$ and,
if $r$ is not `NULL`, $(r, \lceil \mathit{an}/2 \rceil + 1)$ to the
remainder $a - s^2 \le 2s$. Requires $\mathit{an} \ge 1$ and
$a_{an-1} \ne 0$. Inputs shorter than a cutoff depending on the size of
$B$ (43 limbs for $B = 10^{19}$) are converted to binary and use
`flint_mpn_sqrtrem`; longer ones use the Newton square root.
:::

::: function
int radix_sqrt(nn_ptr s, nn_srcptr a, slong an, const radix_t radix)

Sets $(s, \lceil \mathit{an}/2 \rceil)$ to $\lfloor \sqrt a \rfloor$ and
returns 1 if $a$ is a perfect square and 0 otherwise. Requires
$\mathit{an} \ge 1$ and $a_{an-1} \ne 0$.
:::

::: function
void radix_sqrtrem_via_mpn(nn_ptr s, nn_ptr r, nn_srcptr a, slong an, const radix_t radix)
void radix_sqrtrem_newton_karp_markstein(nn_ptr s, nn_ptr r, nn_srcptr a, slong an, const radix_t radix)

The algorithms behind `radix_sqrtrem`, with its interface: conversion to
binary, and the root from a reciprocal square root by one Karp–Markstein
step, with a final correction.
:::

::: function
void radix_rsqrt_1_approx_basecase(nn_ptr res, ulong a, slong n, const radix_t radix)
void radix_rsqrt_1_approx(nn_ptr res, ulong a, slong n, const radix_t radix)

Set $(\mathit{res}, n)$ to the fraction limbs of an approximation of
$1/\sqrt a$ for $2 \le a < B$, with error less than $2 B^{-n}$.
:::

::: function
void radix_rsqrt_approx_basecase(nn_ptr q, nn_srcptr a, slong an, slong n, const radix_t radix)
void radix_rsqrt_approx(nn_ptr q, nn_srcptr a, slong an, slong n, const radix_t radix)

Given $(a, \mathit{an})$ read as a fixed-point number $a \in [B^{-2}, 1)$
with $\mathit{an}$ fraction limbs (one of the top two limbs nonzero), set
$(q, n + 2)$ to an approximation of $1/\sqrt a \in (1, B]$ with $n$
fraction limbs and two integral limbs, with absolute error less than
$4 B^{-n}/\sqrt a$. The range of $a$ is wider than for `radix_inv_approx`
so that inputs with an odd number of limbs can be handled, the scale
factor $B^{1/2}$ not being a limb shift.
:::

::: function
void radix_sqrt_approx_rsqrtmul(nn_ptr q, nn_srcptr a, slong an, slong n, const radix_t radix)
void radix_sqrt_approx(nn_ptr q, nn_srcptr a, slong an, slong n, const radix_t radix)

With $a$ as for `radix_rsqrt_approx`, set $(q, n + 2)$ to an approximation
of $\sqrt a \in [1/B, 1)$ with $n$ fraction limbs and two integral limbs
(normally zero, though the value can round up to 1), with absolute error
less than $4 B^{-n}/\sqrt a$. The error is proportional to $1/\sqrt a$
rather than to $\sqrt a$, because the corrections of the last iteration
are multiplied by an approximation of $1/\sqrt a \le B$.
:::

## Hensel division and $p$-adic square roots {#hensel}

Arithmetic modulo $B^n$, developing results from the least significant
limb upwards. For $b = p$ prime this is arithmetic in $\Z_p$ with
precision $ne$ digits. `mpn_extras` has ports of these algorithms to the
binary radix ([](mpn_extras#hensel)).

::: function
int radix_invmod_bn(nn_ptr res, nn_srcptr x, slong xn, slong n, const radix_t radix)

If $(x, \mathit{xn})$ is invertible modulo $B^n$ (that is, $x_0$ is
coprime to $b$), sets $(\mathit{res}, n)$ to its inverse and returns 1;
otherwise returns 0. Newton iteration on windowed products.
:::

::: function
int radix_divmod_bn(nn_ptr q, nn_ptr rem, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong n, const radix_t radix)

Given $(b, \mathit{bn})$ invertible modulo $B$, sets $(q, n)$ to the
$B$-adic (Hensel) quotient with $qb \equiv a \pmod{B^n}$ and returns 1;
returns 0 if $b$ is not invertible. If `rem` is not `NULL`, the
remainder $(a - qb)/B^n \bmod B^{bn}$ is written to
$(\mathit{rem}, \mathit{bn})$, so that
$a \equiv qb + B^n \mathit{rem} \pmod{B^{n + bn}}$.
:::

::: function
int radix_divmod_bn_1(nn_ptr q, nn_ptr rem, nn_srcptr a, slong an, ulong b, slong n, const radix_t radix)
int radix_divmod_bn_classical(nn_ptr q, nn_ptr rem, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong n, const radix_t radix)
int radix_divmod_bn_karp_markstein(nn_ptr q, nn_ptr rem, nn_srcptr a, slong an, nn_srcptr b, slong bn, slong n, const radix_t radix)

The algorithms behind `radix_divmod_bn`, with its interface: a one-limb
divisor (with `rem` of one limb), the blockwise algorithm (one low and one
high product per block of $\mathit{bn}$ limbs), and the Karp–Markstein
form with a half-precision inverse.
:::

::: function
int radix_sqrtmod_bn(nn_ptr res, nn_srcptr x, slong xn, slong n, const radix_t radix)
int radix_rsqrtmod_bn(nn_ptr res, nn_srcptr x, slong xn, slong n, const radix_t radix)

For a prime digit radix $p$ (with $e \ge 3$ if $p = 2$): if
$(x, \mathit{xn})$ is a unit and a square modulo $p$ (modulo 8 if $p = 2$),
set $(\mathit{res}, n)$ to a square root of $x$ modulo $B^n$, respectively
to the inverse of that square root, and return 1; otherwise return 0. The
square root returned is the one whose residue modulo $p$ is the smallest
nonnegative one.
:::

## Radix conversion {#conversion}

Conversion between radix $B$ and the binary radix $\beta$ in time
$O(M(n) \log n)$ by divide and conquer, with multiplications only (no
divisions). See [](radix/conversion) and
[](algorithms/integers/radix-conversion).

::: function
slong radix_get_mpn(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)
slong radix_get_mpn_basecase(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)
slong radix_get_mpn_divconquer(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)

Convert $(a, \mathit{an})$ from radix $B$ to the binary radix, writing the
limbs to `res` and returning the normalised length of the result. Leading
zero limbs are not counted and may or may not be written. `res` must have
room for the largest integer with $\mathit{an}$ limbs in radix $B$, plus
one extra limb, and must not overlap $a$. The basecase is Horner's rule
(used below 50 limbs).
:::

::: function
slong radix_set_mpn(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)
slong radix_set_mpn_basecase(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)
slong radix_set_mpn_divconquer(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)

Convert the binary $(a, \mathit{an})$ to radix $B$, writing the limbs to
`res` and returning the normalised length of the result. `res` must have
room for `radix_set_mpn_need_alloc(an, radix)` limbs and must not overlap
$a$. The basecase is repeated division by $B$ (used below 40 limbs). The
divide and conquer version caches the most recently used table of powers
(per thread), so that repeated conversions of inputs of the same size, as
when printing many numbers, do not recompute it.
:::

::: function
slong radix_set_mpn_need_alloc(slong n, const radix_t radix)

Returns a number of output limbs that is sufficient for `radix_set_mpn`
with an input of $n$ limbs.
:::

::: type
radix_powers_struct
radix_powers_t

A table of the powers $c^{k_0}, c^{k_1}, \ldots$ used by the divide and
conquer conversions, with $k_0$ the top exponent and
$k_{i+1} = \lceil k_i/2 \rceil$ down to 1, together with their lengths
and their numbers of trailing zero limbs.
:::

::: function
void radix_powers_clear(radix_powers_t powers)

Frees a table of powers.
:::

## String conversion {#strings}

::: function
char * radix_get_str_decimal(char * res, nn_srcptr x, slong n, int negative, const radix_t radix)

Returns the decimal representation of $(x, n)$ as a null-terminated
string, with a leading minus sign if `negative` is 1; for $n = 0$ the
string is `0`. The top limb must be nonzero (if $n > 0$). If `res` is
`NULL`, a new string is allocated with `flint_malloc`; otherwise `res`
must have room for the result and the terminating null. If $b = 10$, this
takes time $O(n)$ (each limb is printed with table lookups for pairs of
digits); otherwise $x$ is converted to binary and then to decimal, in time
$O(M(n) \log n)$.
:::

::: function
char * radix_get_str_sum(char * res, nn_srcptr x, slong n, int negative, int ascending, const radix_t radix)

Returns $(x, n)$ written as a sum over its limbs. For example, 12345678 in
radix $10^3$ is printed as `12 * 10^6 + 345 * 10^3 + 678` with
`negative` and `ascending` both 0, and as
`-(678 + 345 * 10^3 + 12 * 10^6)` with both 1. Memory as for
`radix_get_str_decimal`.
:::

## Integers {#integers}

::: type
radix_integer_struct
radix_integer_t

A signed integer in radix $B$, in a format similar to `mpz_t`: the
fields are `d` (pointer to the limbs), `alloc` (the number of allocated
limbs) and `size`, whose absolute value is the normalised number of limbs
and whose sign is the sign of the integer. The radix is not stored and
must be passed to every function.
:::

::: function
void gr_ctx_init_radix_integer(gr_ctx_t ctx, ulong b, unsigned int exp)

Initialises `ctx` as the generic ring of integers in radix $b^e$ (with
$e = 0$ selecting the largest admissible exponent), whose elements are
`radix_integer_t`. The context owns its `radix_t`.
:::

::: function
void radix_integer_init(radix_integer_t res, const radix_t radix)
void radix_integer_clear(radix_integer_t res, const radix_t radix)
nn_ptr radix_integer_fit_limbs(radix_integer_t res, slong nlimbs, const radix_t radix)

Initialise to zero and free an integer, and make sure that it has room
for `nlimbs` limbs, returning the pointer to its limbs (which the caller
may then fill, setting `size` accordingly).
:::

::: function
int radix_integer_is_normalised(const radix_integer_t x, const radix_t radix)

Returns whether all limbs of $x$ are reduced modulo $B$ and the top limb
is nonzero. Operations on a non-normalised integer have undefined
behaviour.
:::

::: function
void radix_integer_zero(radix_integer_t res, const radix_t radix)
void radix_integer_one(radix_integer_t res, const radix_t radix)
void radix_integer_neg_one(radix_integer_t res, const radix_t radix)
void radix_integer_set(radix_integer_t res, const radix_integer_t x, const radix_t radix)
void radix_integer_set_ui(radix_integer_t res, ulong x, const radix_t radix)
void radix_integer_set_si(radix_integer_t res, slong x, const radix_t radix)
void radix_integer_set_fmpz(radix_integer_t res, const fmpz_t x, const radix_t radix)
void radix_integer_get_fmpz(fmpz_t res, const radix_integer_t x, const radix_t radix)

Assignment and conversion; `radix_integer_set_fmpz` and
`radix_integer_get_fmpz` use the radix conversion functions.
:::

::: function
void radix_integer_rand_limbs(radix_integer_t res, flint_rand_t state, slong max_limbs, const radix_t radix)
void radix_integer_randtest_limbs(radix_integer_t res, flint_rand_t state, slong max_limbs, const radix_t radix)

Set `res` to a random signed integer of at most `max_limbs` limbs, with a
uniform respectively a test distribution of the limbs.
:::

::: function
int radix_integer_is_zero(const radix_integer_t x, const radix_t radix)
int radix_integer_is_one(const radix_integer_t x, const radix_t radix)
int radix_integer_is_neg_one(const radix_integer_t x, const radix_t radix)
int radix_integer_equal(const radix_integer_t x, const radix_integer_t y, const radix_t radix)
int radix_integer_cmp(const radix_integer_t x, const radix_integer_t y, const radix_t radix)
int radix_integer_cmpabs(const radix_integer_t x, const radix_integer_t y, const radix_t radix)
int radix_integer_sgn(const radix_integer_t x, const radix_t radix)

Predicates and comparisons; `cmp` and `cmpabs` return the sign of
$x - y$ and $|x| - |y|$.
:::

::: function
slong radix_integer_size_limbs(const radix_integer_t x, const radix_t radix)
slong radix_integer_size_digits(const radix_integer_t x, const radix_t radix)
slong radix_integer_ssize_limbs(const radix_integer_t x, const radix_t radix)
slong radix_integer_ssize_digits(const radix_integer_t x, const radix_t radix)

The number of limbs, respectively digits, of $|x|$; the `ssize` versions
return it negated if $x$ is negative.
:::

::: function
ulong radix_integer_get_limb(const radix_integer_t x, slong n, const radix_t radix)
ulong radix_integer_get_digit(const radix_integer_t x, slong index, const radix_t radix)

Return the limb, respectively the digit, of $|x|$ at the given
nonnegative position (zero beyond the top).
:::

::: function
void radix_integer_set_limb(radix_integer_t res, const radix_integer_t x, slong index, ulong c, const radix_t radix)
void radix_integer_set_digit(radix_integer_t res, const radix_integer_t x, slong index, ulong c, const radix_t radix)

Set `res` to $x$ with the limb (digit) at the nonnegative position
`index` of $|x|$ replaced by $c < B$ ($c < b$). The sign of $x$ is kept;
if $x$ is zero and $c$ nonzero, the result is positive.
:::

::: function
void radix_integer_neg(radix_integer_t res, const radix_integer_t x, const radix_t radix)
void radix_integer_abs(radix_integer_t res, const radix_integer_t x, const radix_t radix)
void radix_integer_add(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, const radix_t radix)
void radix_integer_sub(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, const radix_t radix)
void radix_integer_mul(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, const radix_t radix)

Negation, absolute value, sum, difference and product.
:::

::: function
void radix_integer_lshift_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)
void radix_integer_lshift_digits(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)
void radix_integer_rshift_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)
void radix_integer_rshift_digits(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)

Set `res` to $x B^n$ ($x b^n$), respectively to $x / B^n$ ($x / b^n$)
truncated towards zero, for $n \ge 0$.
:::

::: function
void radix_integer_addlsh(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, slong v, const radix_t radix)
void radix_integer_sublsh(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, slong v, const radix_t radix)

Set `res` to $x + y b^v$, respectively $x - y b^v$, for any digit shift
$v \ge 0$, without forming the shifted operand: the shift is split into
whole limbs (a pointer offset) and $v \bmod e$ digits (the fused
`radix_addlsh` family). Any aliasing is allowed.
:::

::: function
slong radix_integer_valuation_limbs(const radix_integer_t x, const radix_t radix)
slong radix_integer_valuation_digits(const radix_integer_t x, const radix_t radix)

Return the largest $n$ such that $B^n$ (respectively $b^n$) divides $x$,
that is, the number of trailing zero limbs (digits); 0 for $x = 0$.
:::

::: function
void radix_integer_trunc_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)
void radix_integer_trunc_digits(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)
void radix_integer_mod_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)
void radix_integer_mod_digits(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)
void radix_integer_smod_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)

Reductions modulo $B^n$ (respectively $b^n$): `trunc` reduces $|x|$ and
keeps the sign of $x$; `mod` returns the representative in
$[0, B^n - 1]$ ($[0, b^n - 1]$); `smod` returns the symmetric
representative, in $[-\lfloor B^n/2 \rfloor, \lfloor B^n/2 \rfloor]$ for
odd $B$ and in $[-B^n/2, B^n/2 - 1]$ for even $B$. These take time $O(n)$
plus the time to copy.
:::

::: function
void radix_integer_mullow_limbs(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, slong n, const radix_t radix)

Sets `res` to $xy$ with the absolute value reduced modulo $B^n$, using a
low product.
:::

::: function
int radix_integer_invmod_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)
int radix_integer_divmod_limbs(radix_integer_t res, const radix_integer_t a, const radix_integer_t b, slong n, const radix_t radix)
int radix_integer_rsqrtmod_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)
int radix_integer_sqrtmod_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)

Inverse, quotient $a/b$, inverse square root and square root modulo
$B^n$, by `radix_invmod_bn`, `radix_divmod_bn`, `radix_rsqrtmod_bn` and
`radix_sqrtmod_bn` applied to the absolute values; the sign of the inputs
is carried over to the result. Return 1 on success and 0 if the operand
is not invertible (respectively not a square).
:::

::: function
int radix_integer_div(radix_integer_t q, const radix_integer_t a, const radix_integer_t b, const radix_t radix)
void radix_integer_divexact(radix_integer_t q, const radix_integer_t a, const radix_integer_t b, const radix_t radix)

If $b$ divides $a$, `radix_integer_div` sets $q$ to $a/b$ and returns 1;
otherwise it sets $q$ to zero and returns 0. `radix_integer_divexact`
sets $q$ to $a/b$ assuming that the division is exact (the result is
undefined otherwise).
:::

::: function
void radix_integer_tdiv_qr(radix_integer_t q, radix_integer_t r, const radix_integer_t a, const radix_integer_t b, const radix_t radix)
void radix_integer_fdiv_qr(radix_integer_t q, radix_integer_t r, const radix_integer_t a, const radix_integer_t b, const radix_t radix)
void radix_integer_cdiv_qr(radix_integer_t q, radix_integer_t r, const radix_integer_t a, const radix_integer_t b, const radix_t radix)
void radix_integer_tdiv_q(radix_integer_t q, const radix_integer_t a, const radix_integer_t b, const radix_t radix)
void radix_integer_fdiv_q(radix_integer_t q, const radix_integer_t a, const radix_integer_t b, const radix_t radix)
void radix_integer_cdiv_q(radix_integer_t q, const radix_integer_t a, const radix_integer_t b, const radix_t radix)
void radix_integer_tdiv_r(radix_integer_t r, const radix_integer_t a, const radix_integer_t b, const radix_t radix)
void radix_integer_fdiv_r(radix_integer_t r, const radix_integer_t a, const radix_integer_t b, const radix_t radix)
void radix_integer_cdiv_r(radix_integer_t r, const radix_integer_t a, const radix_integer_t b, const radix_t radix)

Division with remainder $a = qb + r$, where the quotient is rounded
towards zero (`tdiv`), towards $-\infty$ (`fdiv`) or towards $+\infty$
(`cdiv`), as in GMP's `mpz` functions. $b$ must be nonzero.
:::

::: function
int radix_integer_inv(radix_integer_t res, const radix_integer_t x, const radix_t radix)
int radix_integer_is_invertible(const radix_integer_t x, const radix_t radix)

$x$ is invertible in $\Z$ if and only if $x = \pm 1$. `radix_integer_inv`
sets `res` to $1/x = x$ and returns 1 in that case and returns 0
otherwise.
:::

::: function
int radix_integer_sqrt(radix_integer_t s, const radix_integer_t x, const radix_t radix)
void radix_integer_sqrtrem(radix_integer_t s, radix_integer_t r, const radix_integer_t x, const radix_t radix)
int radix_integer_is_square(const radix_integer_t x, const radix_t radix)

`radix_integer_sqrt` sets $s$ to $\lfloor \sqrt x \rfloor$ and returns 1
if $x$ is a perfect square and 0 otherwise; for negative $x$ it returns 0
without modifying $s$. `radix_integer_sqrtrem` sets $s$ and the remainder
$r = x - s^2$, and throws an exception for negative $x$.
`radix_integer_is_square` returns whether $x$ is a perfect square (0 for
negative $x$).
:::

## Division by an invariant word {#div-precomp}

The following inline functions divide words by a fixed word $d \ge 1$ with
a precomputed multiplier, by the method of Granlund and Montgomery
[@GraMon1994] in the variant of Robison [@Rob2005] that needs no
multiplier wider than a word. They are used for the digit shifts of this
module, where $d$ is a power of $b$.

::: type
n_div_precomp_struct
n_div_precomp_t

Precomputed data for division by $d$: the exponent $\ell =
\lfloor \log_2 d \rfloor$ (field `e`), a multiplier `m` and a flag `c`.
If $d$ is a power of two, $m = 0$. Otherwise, with
$q = \lfloor 2^{\ell + \bits}/d \rfloor$, either $m = q + 1$ and $c = 0$
(the *round-up* method, when it is exact for all words) or $m = q$ and
$c = 1$ (the *round-down* method with an increment).
:::

::: type
n_pair_struct

A pair of words (fields `a` and `b`).
:::

::: function
void n_div_precomp_init(n_div_precomp_t pre, ulong d)

Computes the data for division by $d \ge 1$.
:::

::: function
ulong n_div_precomp(ulong x, const n_div_precomp_t pre)
ulong n_div_precomp_bounded(ulong x, const n_div_precomp_t pre)
ulong n_div_precomp_m0(ulong x, const n_div_precomp_t pre)
ulong n_div_precomp_c0(ulong x, const n_div_precomp_t pre)
ulong n_div_precomp_c1(ulong x, const n_div_precomp_t pre)
ulong n_div_precomp_c1_bounded(ulong x, const n_div_precomp_t pre)

Return $\lfloor x/d \rfloor$ for any word $x$: `n_div_precomp` selects
the method from `pre`, and the other functions are the individual methods
for callers that know the case: $x \gg \ell$ for a power of two ($m = 0$),
$\lfloor xm / 2^{\bits} \rfloor \gg \ell$ for $c = 0$, and
$\lfloor (x + 1) m / 2^{\bits} \rfloor \gg \ell$ for $c = 1$, where
$x + 1$ saturates at $\beta - 1$. The `bounded` versions require
$x \ne \beta - 1$ and skip the saturation.
:::

::: function
ulong n_rem_precomp(ulong x, ulong d, const n_div_precomp_t pre)
ulong n_rem_precomp_bounded(ulong x, ulong d, const n_div_precomp_t pre)
ulong n_rem_precomp_m0(ulong x, ulong d, const n_div_precomp_t pre)
ulong n_rem_precomp_c0(ulong x, ulong d, const n_div_precomp_t pre)
ulong n_rem_precomp_c1(ulong x, ulong d, const n_div_precomp_t pre)
ulong n_rem_precomp_c1_bounded(ulong x, ulong d, const n_div_precomp_t pre)
ulong n_divrem_precomp(ulong * r, ulong x, ulong d, const n_div_precomp_t pre)
ulong n_divrem_precomp_bounded(ulong * r, ulong x, ulong d, const n_div_precomp_t pre)
ulong n_divrem_precomp_m0(ulong * r, ulong x, ulong d, const n_div_precomp_t pre)
ulong n_divrem_precomp_c0(ulong * r, ulong x, ulong d, const n_div_precomp_t pre)
ulong n_divrem_precomp_c1(ulong * r, ulong x, ulong d, const n_div_precomp_t pre)
ulong n_divrem_precomp_c1_bounded(ulong * r, ulong x, ulong d, const n_div_precomp_t pre)

The remainder $x \bmod d$, and the quotient with the remainder stored in
`*r`, with the same variants as the quotient functions.
:::

::: function
ulong n_incsat(ulong x)

Returns $x + 1$, or $x$ if $x = \beta - 1$ (saturating increment).
:::
