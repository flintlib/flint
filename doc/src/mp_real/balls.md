---
title: Computing with mp_real balls
summary: How an mp_real_t represents a ball, how to choose the working precision, how radii grow and shrink, and what the arithmetic does not check.
order: 1
---

An `mp_real_t` is a ball $[m \pm r]$ whose midpoint $m$ has a mantissa
of whole limbs and a radix-$B$ exponent ($B = 2^{\bits}$), and whose
radius $r$ is a count of units of the mantissa's last limb. The type is
designed for long computations at high precision in which every
operation should cost as little as possible beyond the limb arithmetic
itself: binary splitting, Newton and AGM iterations, and the evaluation
of constants and elementary functions. Compared with `arb_t`, it has no
infinities or NaN, no special handling of huge exponents, and it rounds
to whole limbs; in exchange an in-place sum or a multiply-add by a word
costs one pass over the limbs of the smaller operand.

## Precision {#precision}

Every operation takes a working precision $n$ in limbs and rounds its
result to about $n$ limbs, or to fewer when an operand is accurate to
fewer. To obtain a relative accuracy of $2^{-p}$ in the basic operations
choose $n$ = `mp_real_prec_bits(p)`, which is
$\lceil p/\bits \rceil + 1$: with one guard limb, sums, products,
quotients and square roots of exact operands have a relative radius of
at most $4 \cdot 2^{-p}$, and of operands accurate to $2^{-p}$ at most
$16 \cdot 2^{-p}$. The elementary functions (`mp_real_exp_bits`, ...)
take the precision in bits directly and choose their limbs internally.

An exact result that fits in $n$ limbs is kept exact, and a longer one
is rounded to $n$ limbs, so that integers computed by exact operations
(as in binary splitting) remain exact until they outgrow the precision.
Inexact products and quotients may keep up to two limbs more than the
accuracy of their operands justifies.

## The radius {#radius}

The radius count `err` is a single limb. It cannot be smaller than one
ulp of the mantissa's last limb unless it is zero, and it cannot reach
$B$ ulps: a radius below one ulp is installed by padding the mantissa
with zero limbs, and a larger one truncates the mantissa. The number of
limbs of an inexact ball is therefore its number of correct limbs, plus
at most one limb of noise. In particular a difference with cancellation
gets shorter, and later operations on it are cheaper.

Rounding errors are rounded up to whole ulps, so every operation can
add one ulp to the radius; and since interval radii accumulate
correlated errors faster than true errors, the arithmetic deliberately
does not keep more than one limb below the accuracy. The functions
`mp_real_add_error_2exp_si` and `mp_real_add_rel_error_2exp_si` add
absolute and relative errors computed by the caller, for example the
tail of a series.

The radius of an operand enters each result through bounds on the
magnitudes of the operands, computed from their leading bits. Divisions
and square roots require operands whose relative radius is below
$2^{-30}$ and throw otherwise: such a ball has hardly any correct bits,
and its quotient would be noise. Convert such balls to `arb_t` with
`mp_real_get_arb` if they must be divided.

## What is not checked {#unchecked}

The arithmetic is rigorous with two exceptions, both for the sake of
speed:

- Exponents are not checked for overflow. A computation whose exponents
  may grow without bound tests its operands with `mp_real_is_safe`
  (exponent at most `MP_REAL_EXP_MAX` limbs in absolute value, $2^{55}$
  on 64-bit machines); one operation on safe operands never overflows.
- Functions with an underscore have preconditions they do not check
  (`_mp_real_submul_bounded` trusts the caller's bound).
  `mp_real_sqrt` and `mp_real_rsqrt` check their operand and throw
  unless it is positive (or, for `mp_real_sqrt`, an exact zero).

## Aliasing and memory {#memory}

All `mp_real_*` functions allow the output to alias any input. In-place
operations are the efficient ones: `mp_real_add(x, x, y, n)` with $y$
inside the window of $x$ costs a limb addition of the length of $y$, and
`mp_real_addmul_ui(x, x, y, c, n)` one `mpn_addmul_1` over $y$. An
`mp_real_t` allocates its mantissa on demand; `mp_real_fit_length`
preallocates. The module keeps process-wide caches of constants and
reduction tables, which `flint_cleanup_master` frees.

## Example {#example}

The program computes $e^{\pi\sqrt{163}}$ at 128 bits, which is famously
close to an integer, and shows that the ball proves that it is not one.
The subtraction of the nearby integer cancels 18 digits, and the result
shrinks to two limbs.

```c example name=mp_real_ramanujan
#include <flint/flint.h>
#include <flint/mp_real.h>
#include <flint/arb.h>

static void show(const char * name, const mp_real_t x)
{
    arb_t a;
    arb_init(a);
    mp_real_get_arb(a, x);
    flint_printf("%s = ", name);
    arb_printn(a, 40, 0);
    flint_printf("   [size %wd, err %wu]\n", x->size, x->err);
    arb_clear(a);
}

int main(void)
{
    mp_real_t pi, r, x;
    slong n = mp_real_prec_bits(128);   /* 3 limbs: 2 + 1 guard limb */

    mp_real_init(pi);
    mp_real_init(r);
    mp_real_init(x);

    mp_real_const_pi4(pi, n, 1);
    mp_real_mul_2exp_si(pi, pi, 2);       /* exact scaling */
    mp_real_set_ui(r, 163);
    mp_real_sqrt(r, r, n);
    mp_real_mul(x, pi, r, n);
    show("pi sqrt(163)", x);

    mp_real_exp_bits(x, x, 128);
    show("exp(pi sqrt(163))", x);

    /* not an integer: the ball contains none */
    flint_printf("unique integer: %d\n", mp_real_unique_integer(r, x));

    /* x - 640320^3 - 744 */
    mp_real_set_ui(r, 640320);
    mp_real_mul(r, r, r, n);
    mp_real_mul_ui(r, r, 640320, n);
    mp_real_sub(x, x, r, n);
    mp_real_set_ui(r, 744);
    mp_real_sub(x, x, r, n);
    show("difference", x);

    mp_real_clear(pi);
    mp_real_clear(r);
    mp_real_clear(x);
    flint_cleanup_master();
    return 0;
}
```

```output
pi sqrt(163) = [40.10916999113251975535008362290414005390 +/- 6.78e-39]   [size 3, err 2]
exp(pi sqrt(163)) = [262537412640768743.99999999999925007260 +/- 4.73e-21]   [size 3, err 525074825281539586]
unique integer: 0
difference = [-7.4992740e-13 +/- 4.73e-21]   [size 2, err 525074825281539586]
```

The radius of the exponential is about $40 \cdot 2^{-128}$ relative to
the result, since the absolute error of the argument becomes a relative
error of $e^x$; this is why the count `err` is large. The exact
products of integers stay exact.

## See also {#see-also}

- [](mp_real/fixed-point) for the fixed-point layer below the balls.
- [](algorithms/floating-point/mp-real-balls) for the representation and
  the error bounds in detail.
- [](arb) for general-purpose ball arithmetic with arbitrary exponents.
