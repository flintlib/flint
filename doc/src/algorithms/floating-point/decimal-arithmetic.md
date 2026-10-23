---
title: Decimal floating-point arithmetic
summary: Correctly rounded decimal arithmetic with limb-aligned exponents in radix 10^19: rounding with sticky bits, sums of far-apart operands, truncated products with interval rounding, division and square roots, decimal radii, conversions to and from binary, and correctly rounded functions by Ziv's strategy with detection of the hard cases.
complexity: Rounding $O(1)$ beyond dropping limbs; sums $O(n)$ for $n$-limb results whatever the exponents; products $O(M(n))$ with rare fallback to a full product; conversions to and from binary $O(M(n)\log n)$.
order: 80
---

The `decimal` module implements decimal floating-point numbers of
arbitrary precision whose operations are correctly rounded to a given
number of significant digits, with unbounded precision and exponents.
Many of its operations and ideas come from IEEE 754-2008 decimal
arithmetic [@IEEE754] and the General Decimal Arithmetic specification
[@Cow2009], but it is not an implementation of either.
The binary floating-point types of FLINT are compared in
[](algorithms/floating-point/overview).

## Representation {#representation}

A nonzero number is $x = (-1)^s M B^v$ where $B = 10^e$ is the limb
radix of the [](radix) module ($e = 19$ on 64-bit machines, the largest
power of ten below $2^{64}$), $M$ is a positive integer stored in
radix-$B$ limbs, and $v$ is an `fmpz` exponent counted in limbs. $M$ is
kept canonical by stripping zero limbs at both ends. The digits of $M$
are thus grouped in limbs of $e$ digits aligned to the exponent, and the
number of significant digits of $x$ is
$$
\operatorname{sd}(x) = \operatorname{digits}(M) - \nu_{10}(M_0),
$$
where $M_0$ is the lowest limb: the trailing zeros of the lowest limb
are an artefact of the alignment and do not count.

Counting the exponent in limbs rather than digits (as `arf` counts bits
and the digit-aligned decimal formats count digits) has three
consequences. Additions align their operands by a whole number of limbs,
with no digit shift. Rounding to $p$ digits drops whole limbs and clears
at most $e - 1$ digits in a single boundary limb. And the decimal
exponent of $x$ is still available in $O(1)$, since only the top limb's
digit count is needed, so the semantics (scientific exponent, ulp,
rounding to $p$ significant digits) do not depend on the alignment. The
price is that a $p$-digit number may need $\lceil p/e \rceil + 1$ limbs,
which matters only for very short numbers. Setting $e = 1$ gives a
digit-aligned representation, used in testing.

## Rounding {#rounding}

For a real $y$ and a precision $p$, let $\operatorname{rnd}_p(y)$ be $y$
rounded to $p$ significant decimal digits in one of the rounding modes
(toward zero or away from zero, toward $\pm\infty$, or to nearest with
ties to even, away from zero or toward zero). Every rounding function
is monotone: $y \le y'$ implies $\operatorname{rnd}_p(y) \le
\operatorname{rnd}_p(y')$. Two simple facts make correct rounding
cheap.

::: lemma Interval rounding {#lem-interval}
If $a \le y \le b$ and $\operatorname{rnd}_p(a) = \operatorname{rnd}_p(b)$,
then $\operatorname{rnd}_p(y) = \operatorname{rnd}_p(a)$.
:::

::: proof
By monotonicity, $\operatorname{rnd}_p(a) \le \operatorname{rnd}_p(y)
\le \operatorname{rnd}_p(b) = \operatorname{rnd}_p(a)$.
:::

::: lemma Sticky bit {#lem-sticky}
Let $y > 0$ have scientific exponent $E$ ($10^E \le y < 10^{E+1}$), and
let $y = Y + \tau$ with $Y$ a multiple of $10^{E - p}$ (one digit beyond
the precision) and $0 \le \tau < 10^{E-p}$. Then $\operatorname{rnd}_p(y)$
is determined by $Y$ and whether $\tau > 0$.
:::

::: proof
The candidates $\lfloor y \rfloor_p$ and $\lceil y \rceil_p$ (the
neighbouring multiples of $10^{E-p+1}$) depend only on $Y$ and on
whether $\tau > 0$, since $y$ is a multiple of $10^{E-p+1}$ exactly when
$\tau = 0$ and the last digit of $Y$ is 0. The directed modes choose
between them by the sign of $y$ alone. The nearest modes compare
$y - \lfloor y \rfloor_p$ with half a unit $5 \cdot 10^{E-p}$: it is
smaller, equal or larger exactly when the last digit $d$ of $Y$
satisfies $d < 5$, ($d = 5$ and $\tau = 0$), or otherwise.
:::

The rounding core `_decfloat_set_round_limbs` takes a mantissa of any
length (with zero limbs allowed at either end) and an optional sticky
flag for an unrepresented tail, and rounds in place: it locates the
rounding position from the digit count of the top limb, drops the limbs
below it, splits the boundary limb by a division by a power of ten (by
precomputed inverse), and propagates an increment if needed. It also
reports whether the result is inexact and bounds the error, which the
ball arithmetic uses as a tight radius.

## Sums {#sums}

The operands of a sum are aligned by their limb exponents and added
exactly, unless the smaller one lies entirely below the rounding
position of the larger. In that case [](#lem-sticky) applies with the
larger operand padded by a limb and the smaller operand reduced to a
signed sticky tail below that limb. The cost is
$O(p)$, independently of the difference of the exponents, so that
$1 + 10^{-10^{12}}$ is rounded as fast as $1 + 10^{-12}$.

## Products {#products}

A product of $n_x$- and $n_y$-limb mantissas is first formed as a high
product: the limbs from position $\ell$ up, computed without the low
partial products, give $r$ with $r \le \lfloor xy/B^\ell \rfloor \le r +
E$ for $E = \min(n_x, n_y, \ell)\,B$. Rounding $r$ and $r + E$ (with the
sticky flag for the discarded low part) gives the correctly rounded
product when the two agree, by [](#lem-interval). They disagree only
when the product lies within $E$ units of a rounding boundary, which
happens with probability about $n/B$ for random operands; the full
product is then computed. The products are those of the [](radix)
module: `radix_mulmid` computes the window of the product directly,
by the classical algorithm for short operands or windows, and otherwise
by the small-prime FFT of `fft_small` (or Kronecker substitution)
applied to the limb sequences as polynomials.

## Division and square root {#division}

A quotient is computed by Euclidean division (`radix_divrem`) of the
numerator, scaled so that the quotient has at least $p + 1$ digits, by
the denominator; the remainder supplies the sticky flag of
[](#lem-sticky), and the quotient is rounded. Square roots use
`radix_sqrtrem` the same way. In the exact ring
$\Z[1/10]$, a quotient exists exactly when the reduced denominator has
only the prime factors 2 and 5, and a square root when the operand is
the square of a decimal fraction.

## Radii {#radii}

A ball radius is a `decmag_t`, $m \cdot 10^k$ with a mantissa $m$ of
exactly `rad_prec` digits ($1 \le$ `rad_prec` $\le 9$ on 64-bit
machines) in one word and an `fmpz` exponent. With $10^{2r+1} < 2^{64}$
for $r$ = `rad_prec`, products of mantissas fit in a word, so each
operation on radii is a few word operations, a division by a power of
ten to renormalise, and an exponent update; the rounding is always in
the direction that preserves the bound. Operands with fewer or more
digits (from a context with a different radius precision) are first
rounded in the safe direction. The ball operations add the propagated
radii and the exact rounding error of the midpoint, rounded up to the
radius precision; for a product,
$$
\operatorname{rad}(xy) = |x_m|\,y_r + |y_m|\,x_r + x_r y_r + |x_m y_m - \operatorname{mid}(xy)|.
$$

## Conversions {#conversions}

Decimal strings convert in linear time, since every limb holds 19
digits. Conversions between binary and decimal use the [](radix)
module, which converts integers between radix $2^{64}$ and radix
$10^{19}$ in $O(M(n)\log n)$ by divide and conquer
([](algorithms/integers/radix-conversion)).

A binary floating-point number $m 2^t$ with $t < 0$ is an exact decimal
number, $m 5^{-t} 10^t$, but with about $0.7|t|$ digits, so correct
rounding to $p$ digits uses exact integer arithmetic on a scaled
truncation (the digits needed plus guard digits and a sticky bit) when
the sizes are moderate. For huge $|t|$ or a huge mantissa the value is
scaled by a power of ten in `arb` arithmetic, $10^k$ being computed in
time essentially linear in the size of $k$, and rounded by Ziv's
strategy; an exactly representable value, which Ziv's strategy could not
decide, is detected first. The conversion from `arb` to a decimal ball
rounds the midpoint correctly and adds the exact conversion error to
the radius. In the other direction, a decimal $M 10^v$ with $v < 0$ is
dyadic only if $5^{-v}$ divides $M$, which is checked cheaply before an
exact conversion to `arf` is attempted.

## Functions {#functions}

The elementary and special functions are evaluated in binary ball
arithmetic by `arb` and `acb` and rounded by Ziv's strategy [@Ziv1991]:
the function is evaluated on the exact argument at increasing working
precision until both endpoints of the resulting ball round to the same
$p$-digit number, which is then the correctly rounded value by
[](#lem-interval). This terminates when the value is not representable,
but no a priori bound on the required precision is known in general
(the table maker's dilemma [@Mul2016]); the implementation gives up with
`GR_UNABLE` at $64p + 1000$ digits, extended by twice the length of the
argument for the elementary functions.

The strategy fails systematically in three situations, which are
detected first, from the sign, the exponent and a few digits of the
argument, in time independent of its length.

- **Exact values.** $\Gamma(n)$, $\log_2(2^n)$, $\sin(\pi k/4)$,
  $\zeta(-2k) = 0$, rational values of the inverse trigonometric
  functions divided by $\pi$, powers with exact results, and so on. For
  example $\sin(\pi x)$ is exact only if $4x$ is an integer, which three
  digits of $x$ decide; a candidate for $\log_2 x$ must have few digits,
  a suitable last digit and a suitable residue modulo a word-size prime
  before $x$ is converted to an integer.
- **Tiny arguments.** If $f(x) = S + t$ where $S$ is an exact decimal
  ($x$, $1$, $\pm 1/x$, $x - 1$, ...) and the tail $t$ has known sign and
  lies more than $p + 4e$ digits below $S$, then $\operatorname{rnd}_p(f(x))$
  is $S$ rounded with $t$ as a signed sticky bit ([](#lem-sticky)). This
  decides $\sin(10^{-10^9})$ in directed rounding modes without any
  evaluation at $10^9$ digits.
- **Large arguments.** Functions converging exponentially to a
  representable limit ($\tanh x \to \pm 1$, $\operatorname{erfc}$,
  $\zeta(s) \to 1$, $\operatorname{expm1}(x) \to -1$) or behaving like
  $1/x$ are handled in the same way, with $S$ the limit.

For complex numbers each component is rounded separately: real and
imaginary arguments are reduced to the real functions (so the exact and
tiny cases carry over), a complex version of the tiny-argument test
determines the sign of each component of $f(z) - S$ from a few Taylor
terms, and the remaining components go through Ziv's strategy with
`acb`.

## Use in FLINT {#use}

- `decfloat_add`, `decfloat_mul`, `decfloat_div`, `decfloat_sqrt` and
  their variants implement the arithmetic, with `_decfloat_set_round_limbs`
  as the rounding core; `_decfloat_add`, `_decfloat_mul` and so on report
  the rounding error for the ball arithmetic.
- `decball_add` and the other ball operations combine them with the
  radius arithmetic of `decmag_t` (`_decmag_add`, `_decmag_mul`, ...).
- `decfloat_set_fmpz_2exp_fmpz`, `decfloat_get_arf`, `decball_set_arb`
  and `decball_get_arb` implement the binary conversions.
- `decfloat_exp` and the other functions, `_decfloat_arb_gr` and
  `_deccfloat_acb_gr` implement Ziv's strategy and the special cases.

::: functions decfloat_mul _decfloat_set_round_limbs decfloat_via_arb _decfloat_round_with_tail :::
