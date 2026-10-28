---
title: GCD, LCM and Bézout coefficients
functions: fmpz_gcd fmpz_gcd3 fmpz_lcm fmpz_xgcd fmpz_xgcd_canonical_bezout fmpz_gcdinv
---

`fmpz_gcd` computes the greatest common divisor of two integers, and
`fmpz_gcd3` and `fmpz_lcm` build on it. For the extended GCD, which also
returns cofactors $a, b$ with $af + bg = \gcd(f, g)$, there are three
functions that differ in how the cofactors are normalised:
`fmpz_xgcd_canonical_bezout` returns the canonical (minimal) cofactors and
is the one to use in new code, `fmpz_xgcd` returns cofactors of another,
older normalisation, and `fmpz_gcdinv` returns only the cofactor of the
smaller operand, reduced modulo the larger, which is what a modular
inverse needs. `fmpz_invmod` ([](fmpz/modular)) is the direct way to
invert modulo an integer.

::: summary algorithms/integers/gcd :::

## Values and signs {#values}

The GCD is always nonnegative, with $\gcd(f, 0) = |f|$ and
$\gcd(0, 0) = 0$; the LCM is nonnegative and 0 if either operand is 0.
All three functions accept any signs and allow aliasing.

## Normalisation of the cofactors {#cofactors}

The cofactors of Bézout's identity are not unique: with $d = \gcd(f, g)$,
$(a + kg/d, b - kf/d)$ is a solution for every $k$.

**`fmpz_xgcd_canonical_bezout`** returns the solution with the smallest
cofactors, which is unique up to the cases below. Writing $d =
\gcd(f, g)$, the result satisfies
$$
|a| < \frac{|g|}{2d}, \qquad |b| < \frac{|f|}{2d},
$$
except in the following cases, where the result is fixed as listed (this
is also the normalisation of GMP's `mpz_gcdext` [@GMP]):

| Case | $(d, a, b)$ |
|---|---|
| $\lvert f \rvert = \lvert g \rvert$ | $(\lvert g \rvert, 0, \operatorname{sgn} g)$ |
| $g = 0$ | $(\lvert f \rvert, \operatorname{sgn} f, 0)$ |
| $f = 0$ | $(\lvert g \rvert, 0, \operatorname{sgn} g)$ |
| $g = \pm 1$ | $(1, 0, g)$ |
| $f = \pm 1$, $g \ne 0, \pm 1$ | $(1, f, 0)$ |
| $\lvert f \rvert = 2d$ | $(d, (d - \lvert g \rvert)/f, \operatorname{sgn} g)$ |
| $\lvert g \rvert = 2d$ | $(d, \operatorname{sgn} f, (d - \lvert f \rvert)/g)$ |

In the last two cases the minimal solution is not unique, since
$|g|/(2d)$ or $|f|/(2d)$ is then an integer or a half-integer that one
cofactor attains.

**`fmpz_xgcd`** normalises the cofactor of the operand of smaller
absolute value instead: if $|f| < |g|$, then $a \operatorname{sgn} f \in
[0, |g|)$ is the cofactor computed by `fmpz_gcdinv` for $|f|$ and $|g|$,
and $b = (d - af)/g$; symmetrically if $|g| < |f|$. For $|f| = |g|$ it
returns $(|f|, \operatorname{sgn} f, 0)$. These cofactors can be larger
than the canonical ones by a factor of up to about $2d$, and the extra
exact division makes the function slower.

**`fmpz_gcdinv`** takes $0 \le f < g$ and returns $d = \gcd(f, g)$ and
$a \in [0, g)$ with $af \equiv d \pmod g$, so that $a = f^{-1} \bmod g$
if $d = 1$. The condition $f < g$ is only checked in builds with
assertions; for $f = 0$ it returns $d = g$ and $a = 0$.

In all three functions the outputs must be distinct objects; the inputs
may alias outputs.

## Algorithms {#algorithms}

Two inline operands, or an inline and a large operand, are handled by
GMP's single-limb GCD `mpn_gcd_1` (a binary GCD, after one division of the
large operand by the small one). Two large operands go to
`flint_mpz_gcd`, which uses a Lehmer-type algorithm with double-limb
steps for small sizes and the subquadratic half-GCD of Möller for large
ones, ported from GMP and built on FLINT's multiplication with
transform reuse in the matrix products; the
extended GCD functions and `fmpz_invmod` use `flint_mpz_gcdext` and
`flint_mpz_invert`, which are built from the same pieces and return the
same cofactors as GMP's `mpz_gcdext` and `mpz_invert`. `fmpz_gcd3` first
computes the GCD of the inline operands, if any, and stops as soon as it
reaches 1; for three large operands it starts with the two shortest.
FLINT's own integer half-GCD, used for continued fractions and rational
reconstruction in `fmpq`, and the partial extended GCD
`fmpz_xgcd_partial` are described with the other algorithms in
[](algorithms/integers/gcd).

## Example {#example}

The cofactors returned by `fmpz_xgcd` and `fmpz_xgcd_canonical_bezout` for
a few pairs, then a GCD and an LCM of Fibonacci numbers (which satisfy
$\gcd(F_m, F_n) = F_{\gcd(m, n)}$) and a modular inverse with
`fmpz_gcdinv`:

```c example name=fmpz_xgcd_variants
#include <flint/flint.h>
#include <flint/fmpz.h>

int main(void)
{
    slong pairs[8][2] = { {240, 46}, {-240, 46}, {12, 18}, {4, 6},
                          {-6, 6}, {0, -5}, {9, -1}, {3, 5} };
    fmpz_t f, g, d, a, b, x;
    int i;

    fmpz_init(f); fmpz_init(g); fmpz_init(d); fmpz_init(a); fmpz_init(b);
    fmpz_init(x);

    flint_printf("   f    g |  xgcd: d   a   b | canonical: d   a   b\n");
    for (i = 0; i < 8; i++)
    {
        fmpz_set_si(f, pairs[i][0]);
        fmpz_set_si(g, pairs[i][1]);
        flint_printf("%4wd %4wd |", pairs[i][0], pairs[i][1]);
        fmpz_xgcd(d, a, b, f, g);
        flint_printf("  %5wd %3wd %3wd |", fmpz_get_si(d), fmpz_get_si(a), fmpz_get_si(b));
        fmpz_xgcd_canonical_bezout(d, a, b, f, g);
        flint_printf("  %5wd %3wd %3wd\n", fmpz_get_si(d), fmpz_get_si(a), fmpz_get_si(b));
    }

    /* gcd, lcm and a modular inverse of large numbers */
    fmpz_fib_ui(f, 300);
    fmpz_fib_ui(g, 200);
    fmpz_gcd(d, f, g);                  /* gcd(F_m, F_n) = F_gcd(m,n) */
    fmpz_lcm(x, f, g);
    flint_printf("gcd(F_300, F_200) = "); fmpz_print(d);
    flint_printf(", lcm has %wu bits\n", fmpz_bits(x));

    fmpz_fib_ui(f, 299);
    fmpz_gcdinv(d, a, g, f);            /* requires 0 <= g < f */
    fmpz_mul(x, a, g);
    fmpz_mod(x, x, f);
    flint_printf("gcd(F_200, F_299) = "); fmpz_print(d);
    flint_printf(", a F_200 mod F_299 = "); fmpz_print(x); flint_printf("\n");

    fmpz_clear(f); fmpz_clear(g); fmpz_clear(d); fmpz_clear(a); fmpz_clear(b);
    fmpz_clear(x);
    return 0;
}
```
```output
   f    g |  xgcd: d   a   b | canonical: d   a   b
 240   46 |      2  -9  47 |      2  -9  47
-240   46 |      2   9  47 |      2   9  47
  12   18 |      6  17 -11 |      6  -1   1
   4    6 |      2   5  -3 |      2  -1   1
  -6    6 |      6  -1   0 |      6   0   1
   0   -5 |      5   0  -1 |      5   0  -1
   9   -1 |      1   0  -1 |      1   0  -1
   3    5 |      1   2  -1 |      1   2  -1
gcd(F_300, F_200) = 354224848179261915075, lcm has 277 bits
gcd(F_200, F_299) = 1, a F_200 mod F_299 = 1
```

For $(12, 18)$ the function `fmpz_xgcd` returns the cofactor $17 \in
[0, 18)$ of the smaller operand, while the canonical solution is
$(-1, 1)$, the special case $|f| = 2d$. The pair $(-6, 6)$ shows the
different conventions for $|f| = |g|$.
