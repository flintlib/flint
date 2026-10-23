---
title: Ball arithmetic with limb-granular radii
summary: The representation of mp_real balls, with the radius counted in ulps of the last limb and the mantissa length following the accuracy; how errors are bounded and installed, and what this costs compared with arb.
complexity: Each operation costs the limb arithmetic of the midpoint plus $O(1)$ word operations for the radius.
order: 30
---

Ball arithmetic represents a real number by a midpoint and a radius and
propagates rigorous bounds through every operation. In `arb`
[@Joh2017] the radius is a `mag_t`, a floating-point number with its own
30-bit mantissa and exponent, and the midpoint an `arf_t` with an
`fmpz` exponent; both are general but cost a few dozen instructions per
operation in bookkeeping. `mp_real` targets long computations at high
precision where operations are counted in limbs anyway (binary
splitting, Newton iterations, evaluation of series), and reduces the
bookkeeping to a few word operations by tying the radius to the
mantissa.

Throughout, $B = 2^{\bits}$.

## Representation {#representation}

::: definition {#def-ball}
An `mp_real_t` with mantissa $(d_0, \ldots, d_{s-1})$, sign $\sigma$,
exponent $e$ and count $c$ represents the ball
$$
\Bigl[(-1)^\sigma\, D\, B^{e - s} \;\pm\; c\, B^{e - s}\Bigr], \qquad
D = \sum_{i<s} d_i B^i,
$$
with $d_{s-1} \ne 0$ when $s > 0$, and $0 \le c < B$. The quantity
$u = B^{e-s}$ is the *ulp* of the ball. For $s = 0$ the midpoint is 0
and the radius is $c\,B^e$.
:::

Two invariants make the representation canonical enough to be cheap:

1. An exact ball ($c = 0$) has $d_0 \ne 0$: exact values carry no low
   zero limbs, so that a small integer stays one limb long however it
   was produced.
2. An inexact ball has $c \ge 1$ and at most one limb of its mantissa
   below the accuracy; the length $s$ is therefore the number of correct
   limbs plus at most one.

::: lemma Relative radius {#lem-rel}
A ball with $s \ge 1$ limbs and count $c$ has relative radius
$c\,u / |m| < c\,B^{1-s}$.
:::

::: proof
$|m| = D\,u \ge d_{s-1} B^{s-1} u \ge B^{s-1} u$.
:::

In particular, at $s = n$ limbs and a count of at most 4, the relative
radius is below $4B^{1-n}$; this is why the working precision
`mp_real_prec_limbs(p)` is $p + 1$ limbs: the guard limb pays for a top
limb that may be as small as 1.

## Installing a bound {#install}

An operation computes its midpoint in a buffer and a bound for its total
error as a count at some anchor (an absolute scale $B^a$). Installing
the bound means choosing the length of the result so that the count
fits in a word.

::: algorithm Installing an error bound {#alg-install}
**Input:** a midpoint $D B^{e-s}$ (exact, $s$ limbs) and a bound
$E \ge 0$ for its error, given as an integer count at scale $B^a$.

**Output:** a ball containing every value within $E B^a$ of the
midpoint, in the representation of [](#def-ball).

1. If $E = 0$, strip the low zero limbs of $D$ and return an exact ball.
2. Express $E B^a$ as a count $c$ of ulps $u = B^{e-s}$, rounding up;
   if the bound is below one ulp, first pad the mantissa with zero limbs
   down to the scale $B^a$, so that $c \ge 1$.
3. If $c \ge B$, let $k \ge 1$ be least with $\lceil c/B^k \rceil + 1 < B$,
   drop the $k$ low limbs of $D$ and set $c \gets \lceil c/B^k \rceil + 1$.
4. Return the ball with the final mantissa and count.
:::

::: lemma {#lem-install}
[](#alg-install) returns a ball that contains the input ball and whose
count is below $B$.
:::

::: proof
Dropping the $k$ low limbs moves the midpoint towards zero by less than
$B^k u$, the new ulp, which the added unit covers; the old count
$c$ is at most $\lceil c/B^k \rceil$ new ulps. Padding with zero limbs
does not change the midpoint. The choice of $k$ makes the count fit.
:::

A consequence is that the radius never has more precision than one
limb: a count of $c$ at a given ulp is accurate to $\lceil \log_2 c
\rceil$ bits at most. This is enough, because a radius is an upper bound
and only its order of magnitude matters, and it saves the second
exponent that `mag_t` carries. Rounding a bound up to whole ulps costs
up to one ulp per operation, which a guard limb absorbs.

## Bounds in limb arithmetic {#bounds}

Internally, bounds are two-limb counts at a limb anchor: the product of
a count and a magnitude is one `umul_ppmm`, and the sum of two counts
one `add_ssaaaa` after aligning the lower anchor to the higher by a
division by a power of $B$, rounding up. The magnitudes of operands are
bounded by their leading $\bits$ bits, read across a limb boundary;
bounding through the top limb alone would overstate a value like
$1 + \varepsilon$, whose top limb is 1, by a factor 2, and with radii
held as whole counts such pessimism turns into truncated limbs.
Divisions and square roots bound their error through doubles (the
derivative of the function at the leading bits of the operand) and
convert the result to a count when it is installed; the doubles are kept
in their normal range by explicit scaling, so that no floating-point
exception is raised.

## Sums {#sums}

A sum is formed in the destination buffer. A difference is taken from
the operand of larger magnitude, decided up front from the exponents
and top limbs, so that no negation pass follows. The bottom limb of an
inexact operand is its noise floor: no ball holds information below it,
so a sum takes nothing below the floor of an inexact operand, and the
discarded tail of the other operand costs one unit there.

::: proposition In-place sums {#prop-inplace}
Let $x$ be inexact and $y$ lie within the window of $x$ (its exponent at
most that of $x$). Then $x := x \pm y$ costs $O(s_y)$ limb operations
plus $O(1)$, independently of the length of $x$, unless the result
changes the length of $x$.
:::

::: proof
The limbs of $y$ above the floor of $x$ are added to (or subtracted
from) the corresponding limbs of $x$; the carry or borrow propagates
into the higher limbs of $x$ only while it is nonzero, which costs
amortised $O(1)$. The part of $y$ below the floor is replaced by one
unit added to the count. Normalisation (a new top limb or a vanished
one) and the installation of the count are constant time unless the
top limb changes, in which case the mantissa is shifted by one limb.
:::

This is the property that makes rectangular splitting and
multiply-accumulate loops efficient: `mp_real_addmul_ui` in place costs
one `mpn_addmul_1` over $y$.

## Products and quotients {#products}

For a product $ab$ at precision $n$, both operands are first truncated
to $n + 2$ limbs; the dropped tail of $a$, below one unit of its lowest
kept limb, enters the bound multiplied by $|b|$ and vice versa. A long
exact operand therefore costs no more than the precision asks for. The
radius of the product is

$$
|a|\,\rho_b + |b|\,\rho_a + \rho_a\rho_b + (\text{truncations}) +
(\text{rounding of the product}),
$$

with $|a|$, $|b|$ bounded through their leading bits. Short products are
formed in full; longer ones as middle products (`flint_mpn_mulmid`),
which compute only the window of limbs that reaches the result, with a
bounded error for the discarded low part
([](algorithms/integers/short-products)).

Quotients and roots use the fixed-point Newton iterations
(`_mp_real_div_newton`, `_mp_real_rsqrt_newton`, `_mp_real_sqrt_newton`,
see [](algorithms/integers/division)), whose errors are bounded
relative to the result. The operand radii enter through the derivative:
for $a/b$, $\rho_a/|b| + |a|\rho_b/|b|^2$, multiplied by $1.01$ to cover
the correction $|b| - \rho_b$ in the denominator. That factor is valid
because the relative radius of $b$ is required to be below $2^{-30}$;
`mp_real_div` and `mp_real_sqrt` check this and throw otherwise. A ball
that wide has no correct mantissa limbs, so a quotient by it is a usage
error rather than something worth bounding precisely.

## Integers from balls {#unique-integer}

`mp_real_unique_integer` extracts the integer contained in a ball, as
needed after an evaluation of the partition function, whose result is
known to be an integer.

::: lemma {#lem-unique}
Let $m = I + f$ with $I \in \Z$, $0 \le f < 1$, and $0 \le \rho < 1$.
The ball $[m \pm \rho]$ contains exactly one integer if and only if
$[f + \rho \ge 1] = [f - \rho > 0]$, and that integer is then
$I + [f + \rho \ge 1]$. (Here $[\cdot]$ is 1 for a true statement and
0 otherwise.)
:::

::: proof
The integers in the ball are those from $\lceil m - \rho \rceil$ to
$\lfloor m + \rho \rfloor$. Since $f + \rho \in [0, 2)$ and
$f - \rho \in (-1, 1)$, $\lfloor m + \rho \rfloor = I + [f + \rho \ge 1]$
and $\lceil m - \rho \rceil = I + [f - \rho > 0]$. There is exactly one
integer when the two are equal.
:::

For a positive midpoint whose fraction has $k$ limbs and a radius below
one ulp of the integer part, the two conditions are the carry out of
$f + \rho$ and the borrow (or zero result) of $f - \rho$ in $k$-limb
arithmetic. Both are decided by scanning $f$ from its top limb until a
limb that is neither all zeros nor all ones, which is almost always the
first one; a negative midpoint is handled by symmetry. The integer part
is then moved down in place (and incremented if needed), without the
size limits of `fmpz_t`.

## Use in FLINT {#use}

All `mp_real_*` ball functions use this representation:
`mp_real_add`, `mp_real_mul`, `mp_real_div`, `mp_real_sqrt` and
their variants implement the rules above. The binary splittings of
the module (`mp_real_hypgeom_series`, the bit-burst steps) rely on exact
balls staying exact and short and on inexact ones losing limbs as they
lose accuracy; the elementary functions use balls around their
fixed-point kernels so that every error is accounted for
([](algorithms/floating-point/elementary-functions)).

::: functions mp_real_add mp_real_mul mp_real_div mp_real_unique_integer :::
