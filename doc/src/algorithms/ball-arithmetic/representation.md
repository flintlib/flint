---
title: Representation of radii and midpoints
summary: The formats of mag_t (30-bit mantissa, upper-bound rounding) and arf_t (variable-length mantissa, correct rounding), with proofs that the basic mag operations give upper bounds, and how a rounding error becomes part of the radius.
complexity: A mag_t operation costs a few word operations; an arf_t operation on n-limb mantissas costs O(n) for addition and M(n) for multiplication, plus O(1) word operations for rounding.
order: 10
---

A ball $[m \pm r]$ in FLINT has an `arf_t` midpoint and a `mag_t`
radius. The two types are designed for different jobs. The midpoint
carries the value and needs arbitrary precision with correct rounding;
the radius carries an error bound, for which a few significant bits
suffice as long as every operation rounds upwards. The design is
discussed in [@Joh2017].

## The mag_t format {#mag}

A finite nonzero `mag_t` stores a 30-bit integer mantissa $a$ with
$2^{29} \le a < 2^{30}$ and an `fmpz` exponent $e$, and represents

$$
x = a \cdot 2^{e - 30} \in [2^{e-1}, 2^e).
$$

Zero is $a = 0$, $e = 0$ and $+\infty$ is $a = 0$, $e \ne 0$. There are
no negative numbers and no NaN. The exponent is an `fmpz`, so it is a
plain `slong` stored inline unless it exceeds about $2^{62}$ in absolute
value; all functions handle bignum exponents, but the common case costs
no memory allocation and a comparison of exponents is a comparison of
two words.

Thirty bits are chosen so that a product of two mantissas fits in one
64-bit word, with room for the carry of a few additions: the product of
two mantissas is computed as `MAG_FIXMUL(a, b)` $= \lfloor ab / 2^{30}
\rfloor$ with a single machine multiplication, also on 32-bit systems
(through a 64-bit type). Error bounds are needed only to a few digits, so
30 bits lose nothing in practice, and they make `mag_t` arithmetic many
times faster than `arf_t` arithmetic at any precision.

Every operation rounds in a fixed direction: the functions without
suffix return upper bounds and those with suffix `_lower` lower bounds.
Rather than determining whether the exact result is representable, the
upper-bound functions simply add one unit in the last place. The results
are not optimal, but they are always on the safe side, and the code has
no branches for exactness.

::: algorithm Multiplication {#alg-mag-mul}
**Input:** $x = a 2^{e-30}$, $y = b 2^{f-30}$ with $2^{29} \le a, b <
2^{30}$.

**Output:** an upper bound $z \ge xy$ in `mag_t` format.

1. $c \gets \lfloor ab / 2^{30} \rfloor + 1$, $g \gets e + f$.
2. If $c < 2^{29}$, set $c \gets 2c$ and $g \gets g - 1$.
3. Return $z = c \cdot 2^{g - 30}$.
:::

::: proposition {#prop-mag-mul}
[](#alg-mag-mul) returns a normalised `mag_t` with $xy < z \le xy +
2^{g - 29}$, that is, with an error of at most two units in the last
place.
:::

::: proof
We have $xy = ab \cdot 2^{e + f - 60}$, and $c_1 = \lfloor ab/2^{30}
\rfloor + 1 > ab/2^{30}$, so $c_1 2^{e+f-30} > xy$. Doubling $c_1$ and
decreasing the exponent does not change the value. Since $2^{58} \le ab
\le (2^{30} - 1)^2 < 2^{60} - 2^{31}$, we get $2^{28} \le
\lfloor ab/2^{30} \rfloor \le 2^{30} - 2$, so $2^{28} < c_1 < 2^{30}$,
and after step 2 the mantissa lies in $[2^{29}, 2^{30})$. The excess
$c_1 - ab/2^{30}$ is at most 1 unit of $2^{e+f-30}$, which is at most 2
units of the final exponent.
:::

Addition of $x = a 2^{e-30}$ and $y = b 2^{f-30}$ with $e \ge f$ adds $a$
and $\lfloor b / 2^{e-f} \rfloor + 1$ (or just 1 if $e - f \ge 30$, since
then $y < 2^{e-30}$). The sum is at most 31 bits; if it overflows the
mantissa, it is shifted right by one bit and the lost bit is added back,
which again rounds up. The lower-bound versions omit the added units and
truncate. Division computes $\lfloor a 2^{30} / b \rfloor + 1$ with one
word division. Square roots, logarithms and exponentials are computed with
`double` arithmetic on the mantissa and multiplied by a correction factor
such as $1 \pm 10^{-13}$ before conversion back with an extra unit in the
last place; this is valid under the weak assumptions on floating-point
arithmetic listed on [](guide/types/ball-arithmetic#conventions), because the
relative error of each `double` operation is far below the correction
factor. The logarithm is the only transcendental function needed for
`mag_t` bounds, and since the C standard says nothing about the accuracy
of `log`, FLINT uses its own implementation with a proven error bound
(`mag_d_log_upper_bound`).

The `mag_fast_` functions assume that all exponents are small and
finite, skip the checks for special values and bignum exponents, and
are used in inner loops such as `arb_dot` and `arb_mul`.

## The arf_t format {#arf}

A finite nonzero `arf_t` stores a sign, an `fmpz` exponent $e$ and a
mantissa of $n$ limbs $d_{n-1}, \ldots, d_0$ with the top bit of
$d_{n-1}$ set and $d_0 \ne 0$; it represents

$$
x = \pm 0.d_{n-1} d_{n-2} \ldots d_0 \cdot 2^e \in \pm [2^{e-1}, 2^e).
$$

Mantissas of one or two limbs are stored inside the structure, longer
ones in a heap-allocated array that grows as needed. The normalisation
$d_0 \ne 0$ means that trailing zero limbs are never stored: a number has
a unique representation, an integer such as 3 has a one-limb mantissa
regardless of the precision of the computation that produced it, and
`arf_bits` is cheap. Zero, $\pm\infty$ and NaN have a size field of zero
and are distinguished by special exponent values.

The precision is not a property of an `arf_t`. Each operation takes a
precision $p$ and rounds its result to at most $p$ bits, correctly, in a
given direction; the result can be shorter. This makes it possible to mix
exact and approximate values freely, which ball arithmetic needs
constantly: an exact integer, a power of two or a short dyadic number
stays cheap however high the working precision is.

Rounding a mantissa of $k$ bits to $p < k$ bits (`_arf_set_round_mpn`)
keeps the top $p$ bits and inspects the remaining ones only to decide
whether they are all zero (for the inexact flag) and, for rounding to
nearest, whether they are above, below or equal to half a unit. Rounding
up can carry into a new top bit, in which case the mantissa becomes a
power of two and the exponent increases. With `ARF_RND_DOWN`, the mode
used for ball midpoints, no carry can occur, so this is the cheapest mode.

Addition aligns the mantissas and adds them, but it never forms an
exact sum longer than about $p$ bits plus a few limbs: if the smaller
operand lies entirely below the last bit of the result, it only
influences the rounding, and a sticky bit suffices. Multiplication
computes the full product of the mantissas, except at intermediate sizes
where `arf_mul_rnd_down` uses the short product of MPFR [@BZ2011].
Division and square root use a quotient or root with 32 extra bits, check
exactness by multiplication when the extra bits are zero, and, in the
relaxed modes `ARF_RND_FAST` and `ARF_RND_ACCURATE` at high precision,
Newton iteration from [mp_real](mp_real).

## From rounding error to radius {#ulp}

Let $\tilde m$ be the midpoint of a result rounded to $p$ bits with
`ARF_RND_DOWN`, and let $2^{e-1} \le |\tilde m| < 2^e$. Truncation to $p$
bits changes the exact value $m$ by less than one unit in the last place
of $\tilde m$, so

$$
|m - \tilde m| < 2^{e - p}.
$$

`arf_mag_add_ulp(r, r, m, p)` adds $2^{e - p}$ to the radius $r$ (with
`mag_add_2exp_fmpz`, which rounds upward). This is done only if the
operation reported an inexact result, so that operations on exact input
with an exact result give exact balls. The relaxed modes of division and
square root have an error of at most $0.5625$ units of $2^{e' - p}$,
where $e'$ is the exponent of the exact result, and return either a value
with the same exponent or exactly $2^{e'}$, so the same bound applies.
They do not report exactness reliably, however; `arb_div` and the square
root functions rely on the flag anyway, which is a defect described on
`arb_div`.

## Use in FLINT {#use}

- `mag_t` ([mag](mag)) is the radius of `arb_t`, `acb_t` and of the
  coefficients of `arb_poly`, `arb_mat` and so on, and the type of almost
  all error bounds computed inside algorithms.
- `arf_t` ([arf](arf)) is the midpoint type, and is used directly where
  correctly rounded floating-point arithmetic is wanted, for example in
  `arb_fpwrap` and in approximate algorithms with `arf` or `acf` numbers.
- `arf_mag_add_ulp`, `arf_mag_set_ulp` and `arf_mag_fast_add_ulp` turn
  the inexact flag of a midpoint operation into a radius contribution.

::: functions mag_mul mag_add arf_add arf_mul arf_mag_add_ulp :::
