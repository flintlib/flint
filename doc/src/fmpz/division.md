---
title: Division and rounding conventions
functions: fmpz_tdiv_qr fmpz_ndiv_qr fmpz_tdiv_q fmpz_fdiv_r fmpz_tdiv_q_2exp fmpz_tdiv_ui fmpz_divexact fmpz_divexact2_uiui fmpz_div
---

For integers $g$ and $h \ne 0$ there are many pairs $(q, r)$ with
$g = qh + r$ and $|r| < |h|$; a division function picks one by rounding
the rational quotient $g/h$ to an integer. The `fmpz` functions name the
rounding in their prefix, following GMP: `tdiv` truncates (rounds towards
zero), `fdiv` takes the floor, `cdiv` the ceiling, and `ndiv` rounds to
the nearest integer. A suffix `_qr` returns quotient and remainder, `_q`
the quotient only and `_r` the remainder only; `_ui` and `_si` take a
word-size divisor and `_2exp` divides by a power of two. Exact division,
when the remainder is known to be zero, is a separate and faster
operation (`fmpz_divexact`), and `fmpz_div` checks exactness.

::: summary algorithms/integers/division :::

## The four conventions {#conventions}

| Prefix | Quotient $q$ | Sign of the remainder $r = g - qh$ | Range of $r$ |
|---|---|---|---|
| `tdiv` | $\operatorname{trunc}(g/h)$ | that of $g$ | $\lvert r \rvert < \lvert h \rvert$ |
| `fdiv` | $\lfloor g/h \rfloor$ | that of $h$ | $0 \le r < h$ or $h < r \le 0$ |
| `cdiv` | $\lceil g/h \rceil$ | opposite to that of $h$ | $-h < r \le 0$ or $0 \le r < -h$ |
| `ndiv` | nearest integer to $g/h$, ties to even | either | $\lvert r \rvert \le \lvert h \rvert / 2$ |

A zero remainder satisfies every sign condition. Which to choose:

- `fdiv` with $h > 0$ is reduction modulo $h$: the remainder lies in
  $[0, h)$, as returned by `fmpz_mod` (which also accepts negative $h$,
  using $|h|$). `fmpz_fdiv_q_2exp` is the arithmetic right shift and
  `fmpz_fdiv_r_2exp` the reduction modulo $2^e$.
- `tdiv` matches C's `/` and `%` and is symmetric under $g \mapsto -g$.
- `cdiv` gives the smallest $q$ with $qh \ge g$ for $h > 0$, as needed for
  rounding up sizes or bounds.
- `ndiv` gives a remainder of minimal absolute value, as needed for
  balanced representations; `fmpz_smod` gives the corresponding symmetric
  residue, but with ties resolved to $+|h|/2$ instead of by the parity of
  the quotient.

The functions `fmpz_tdiv_ui`, `fmpz_fdiv_ui` and `fmpz_cdiv_ui` return
the absolute value of the remainder as an `ulong`; only `fmpz_fdiv_ui`
(for which it is the remainder itself) is a residue in $[0, h)$.

The quotients by $2^e$ are shifts of the two's complement representation
for `fdiv` and of the absolute value for `tdiv`, so `fmpz_tdiv_q_2exp`
and `fmpz_fdiv_q_2exp` differ for negative $g$ that are not multiples of
$2^e$.

Division by zero raises an exception, except in `fmpz_mod` and
`fmpz_smod`, where it is an unchecked precondition, and in `fmpz_div`,
which reports it as an inexact division unless $g = 0$.

## Exact division {#exact}

`fmpz_divexact` divides $g$ by a divisor $h$ that is known to divide it.
The quotient is then determined by the low limbs of $g$ and $h$ alone, and
FLINT computes it by 2-adic (Hensel) division from the least significant
limb upwards, avoiding the correction steps of Euclidean division. If
$h$ does not divide $g$ the result is meaningless; nothing is checked.
Typical uses are dividing out a known content or GCD, the last step of an LCM, binomial
coefficients, and fraction-free elimination.

`fmpz_div` (also available under its older name `fmpz_divides`) is the
checked version: it returns 1 and the quotient if $h$ divides $g$, and 0
otherwise. `fmpz_divisible` only answers the question, which for large
operands is decided from the low limbs by a 2-adic method without forming
the whole quotient.

## Algorithms {#algorithms}

Inline operands are divided with the C operators and a sign correction.
A large dividend with a divisor of one word is divided by one or two
hardware divisions if it has at most two limbs, and otherwise by GMP's
`mpz_tdiv_q_ui` and its relatives, built on GMP's `mpn_divrem_1`
(division by a single limb with a precomputed inverse,
[](algorithms/arithmetic/moller-granlund)). Two large operands
are divided by `_flint_mpn_tdiv_qr`, FLINT's replacement of GMP's
division. Depending on the sizes it uses schoolbook division with a 3-by-2 inverse, the
divide-and-conquer division of Burnikel and Ziegler, block division with
a shared approximate inverse for long quotients, or Newton iteration
(Karp–Markstein division) once divisor and quotient have about 950 limbs.
Exact division (`fmpz_divexact`, `flint_mpn_divexact`) uses the 2-adic
analogues of the same algorithms. [](algorithms/integers/division)
describes them; `fmpz_fdiv_qr_preinvn` exposes the precomputed-inverse
division for repeated divisions by the same divisor.

## Example {#example}

All four conventions for $\pm 7$ divided by $\pm 2$ (each pair is
$q$, $r$), shifts of a negative number, and exact and checked division
of $30!$:

```c example name=fmpz_division_conventions
#include <flint/flint.h>
#include <flint/fmpz.h>

int main(void)
{
    slong g[4] = { 7, -7, 7, -7 }, h[4] = { 2, 2, -2, -2 };
    fmpz_t a, b, q, r;
    int i;

    fmpz_init(a); fmpz_init(b); fmpz_init(q); fmpz_init(r);

    flint_printf("  g   h |  tdiv   |  fdiv   |  cdiv   |  ndiv\n");
    for (i = 0; i < 4; i++)
    {
        fmpz_set_si(a, g[i]);
        fmpz_set_si(b, h[i]);
        flint_printf("%3wd %3wd |", g[i], h[i]);
        fmpz_tdiv_qr(q, r, a, b);
        flint_printf(" %2wd %2wd  |", fmpz_get_si(q), fmpz_get_si(r));
        fmpz_fdiv_qr(q, r, a, b);
        flint_printf(" %2wd %2wd  |", fmpz_get_si(q), fmpz_get_si(r));
        fmpz_cdiv_qr(q, r, a, b);
        flint_printf(" %2wd %2wd  |", fmpz_get_si(q), fmpz_get_si(r));
        fmpz_ndiv_qr(q, r, a, b);
        flint_printf(" %2wd %2wd\n", fmpz_get_si(q), fmpz_get_si(r));
    }

    /* shifts: division by 2^3 of -100 */
    fmpz_set_si(a, -100);
    fmpz_tdiv_q_2exp(q, a, 3);
    fmpz_fdiv_r_2exp(r, a, 3);
    flint_printf("tdiv_q_2exp(-100, 3) = %wd, fdiv_r_2exp(-100, 3) = %wd\n",
                 fmpz_get_si(q), fmpz_get_si(r));
    fmpz_fdiv_q_2exp(q, a, 3);
    flint_printf("fdiv_q_2exp(-100, 3) = %wd, tdiv_ui(-100, 8) = %wu\n",
                 fmpz_get_si(q), fmpz_tdiv_ui(a, 8));

    /* exact and checked division of 30! by 2^26 3^14 */
    fmpz_fac_ui(a, 30);
    fmpz_one_2exp(b, 26);
    fmpz_mul_ui(b, b, 4782969);         /* 3^14 */
    fmpz_divexact(q, a, b);
    flint_printf("30!/(2^26 3^14) = "); fmpz_print(q); flint_printf("\n");
    fmpz_mul_ui(b, b, 3);
    flint_printf("divisible by 3^15: %d\n", fmpz_div(q, a, b));

    fmpz_clear(a); fmpz_clear(b); fmpz_clear(q); fmpz_clear(r);
    return 0;
}
```
```output
  g   h |  tdiv   |  fdiv   |  cdiv   |  ndiv
  7   2 |  3  1  |  3  1  |  4 -1  |  4 -1
 -7   2 | -3 -1  | -4  1  | -3 -1  | -4  1
  7  -2 | -3  1  | -4 -1  | -3  1  | -4 -1
 -7  -2 |  3 -1  |  3 -1  |  4  1  |  4  1
tdiv_q_2exp(-100, 3) = -12, fdiv_r_2exp(-100, 3) = 4
fdiv_q_2exp(-100, 3) = -13, tdiv_ui(-100, 8) = 4
30!/(2^26 3^14) = 826385373016328125
divisible by 3^15: 0
```

For $7/2 = 3.5$ the nearest-integer rounding has a tie and chooses the
even quotient 4; for $-7/2$ it chooses $-4$. The exponents 26 and 14 are
the 2-adic and 3-adic valuations of $30!$, so the exact quotient is odd
and not divisible by 3.
