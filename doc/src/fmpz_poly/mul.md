---
title: Multiplying polynomials
functions: fmpz_poly_mul fmpz_poly_mullow fmpz_poly_mulmid fmpz_poly_mulhigh_n _fmpz_poly_mulhigh fmpz_poly_sqr fmpz_poly_sqrlow
---

`fmpz_poly_mul` multiplies two polynomials in $\Z[x]$,
`fmpz_poly_mullow` computes the product modulo $x^n$,
`fmpz_poly_mulmid` an arbitrary window of coefficients, and
`fmpz_poly_mulhigh_n` the top coefficients. `fmpz_poly_sqr` and
`fmpz_poly_sqrlow` are the corresponding squares. All of them call one
dispatcher, `_fmpz_poly_mulmid`, which chooses among a dozen algorithms
from the lengths, the window and the sizes of the coefficients, so
there is rarely a reason to call an algorithm directly.

::: summary algorithms/polynomials/multiplication-integers :::

## Choosing the algorithm {#selection}

Polynomials over $\Z$ vary in two dimensions, the length $n$ and the
coefficient size $b$, and the cost of a product depends on both: it is
roughly that of multiplying integers of $nb$ bits when $n$ is large,
and of $n^2$ (or fewer) products of $b$-bit integers when $n$ is small.
Writing $b_1, b_2$ for the bit sizes of the coefficients of the two
operands, $s = b_1 + b_2$, and $n$ for the length of the shorter
operand (after discarding coefficients that do not contribute to the
window), `_fmpz_poly_mulmid` proceeds as follows.

1. A constant operand is handled as a scalar multiplication.
2. For middle windows (touching neither end of the product) with $n
   \ge 16$ and $s \approx 2000$–$4400$ bits or $s \ge 120\,000$ bits,
   a Schönhage–Strassen cyclic convolution of at most three quarters of
   the truncated length is used if one exists, since its wrap-around
   then falls outside the window.
3. With `fft_small`, products with $n$ from about 80 to 130 on (the
   cutoff depends on $s$ and the window) are computed by the
   multimodular small-prime FFT `_fmpz_poly_mul_mid_default_mpn_ctx`,
   which works directly on the coefficients when $s$ plus the bits of
   $n$ is below about 397, and declines otherwise.
4. Short products with coefficients of at most 62 bits whose result
   fits in one or two words are accumulated directly in machine words.
5. Very short operands use the classical algorithm, and operands whose
   coefficient sizes differ by a factor of at least 8 use the classical
   algorithm or the packed-limb functions `_fmpz_poly_mulmid_mpn_bits`,
   which only pay for the actual sizes.
6. For $s < 3000$ bits: the classical algorithm up to a few terms
   (between 2 and 15, depending on $s$ and the window), then the
   packed-limb classical and Karatsuba algorithms up to 20–80 terms.
7. For $s \ge 3000$ bits: Karatsuba with Toom–Cook base cases
   (`_fmpz_poly_mulmid_toom_karatsuba`) up to 4–18 terms; for truncated
   products with $s \ge 16\,000$ bits, classical multiplication in the
   transform space of `fft_small` up to 21–87 terms.
8. Otherwise [Kronecker substitution](algorithms/polynomials/kronecker)
   (`_fmpz_poly_mulmid_KS`), whose single integer product is computed
   by `flint_mpn_mul` and therefore by `fft_small` when it is large;
   Schönhage–Strassen (`_fmpz_poly_mulmid_SS`) is used instead for
   squares of length 8 to 75 with $s$ between 800 and 4000 bits, for
   huge products ($s \ge 80\,000$ bits and $(\mathit{len1} +
   \mathit{len2}) s > 4 \cdot 10^9$), and with four or more threads for
   long products with large coefficients. Builds without `fft_small`
   choose between KS and SS by the number of limbs per coefficient.

The boundaries were measured on x86-64 for random coefficients of
uniform size. Coefficients of very different sizes within one operand
are detected for short operands: the bounds are then reduced to the
average size, which favours the direct algorithms. Squares are
detected when both operands are the same array with the same length,
and are up to twice as fast.

The truncated and middle products are not just truncations of the full
product. The classical and packed-limb algorithms only compute the
products that contribute to the window; Karatsuba uses Mulders' short
product for low windows and the middle product of Hanrot, Quercia and
Zimmermann for middle windows; Kronecker substitution computes only the
needed limbs of the integer product with `flint_mpn_mulmid`; and the
FFTs use truncated or cyclic transforms. A middle product of a $2m$ by
$m$ product costs about as much as an $m$ by $m$ product
([](algorithms/polynomials/middle-product)), which is what Newton
iteration needs.

## Examples {#examples}

The parts of one product. The low product ends in a zero coefficient
and is normalised; the high product has the full length, and only its
top 3 coefficients are meaningful.

```c example name=fmpz_poly_mul_parts
#include <flint/flint.h>
#include <flint/fmpz_poly.h>

int main(void)
{
    fmpz_poly_t a, b, c;

    fmpz_poly_init(a);
    fmpz_poly_init(b);
    fmpz_poly_init(c);

    fmpz_poly_set_str(a, "4  1 -2 3 4");       /* 4x^3 + 3x^2 - 2x + 1 */
    fmpz_poly_set_str(b, "3  5 0 -7");         /* -7x^2 + 5 */

    fmpz_poly_mul(c, a, b);
    flint_printf("a*b          = "); fmpz_poly_print_pretty(c, "x"); flint_printf("\n");
    fmpz_poly_mullow(c, a, b, 4);
    flint_printf("mullow 4     = "); fmpz_poly_print_pretty(c, "x"); flint_printf("\n");
    fmpz_poly_mulmid(c, a, b, 2, 5);
    flint_printf("mulmid [2,5) = "); fmpz_poly_print_pretty(c, "x"); flint_printf("\n");
    fmpz_poly_mulhigh_n(c, a, b, 3);
    flint_printf("mulhigh_n 3  = top coefficients %wd %wd %wd of length %wd\n",
        fmpz_poly_get_coeff_si(c, 5), fmpz_poly_get_coeff_si(c, 4),
        fmpz_poly_get_coeff_si(c, 3), fmpz_poly_length(c));
    fmpz_poly_sqr(c, b);
    flint_printf("b^2          = "); fmpz_poly_print_pretty(c, "x"); flint_printf("\n");

    fmpz_poly_clear(a);
    fmpz_poly_clear(b);
    fmpz_poly_clear(c);
    return 0;
}
```
```output
a*b          = -28*x^5-21*x^4+34*x^3+8*x^2-10*x+5
mullow 4     = 34*x^3+8*x^2-10*x+5
mulmid [2,5) = -21*x^2+34*x+8
mulhigh_n 3  = top coefficients -28 -21 34 of length 6
b^2          = 49*x^4-70*x^2+25
```

All algorithms compute the same product. The program multiplies random
polynomials of length 60 with 1000-bit coefficients with each of them,
and checks a truncated and a middle product against the full product:

```c example name=fmpz_poly_mul_algorithms
#include <flint/flint.h>
#include <flint/fmpz_poly.h>

int main(void)
{
    flint_rand_t state;
    fmpz_poly_t a, b, c, d;

    flint_rand_init(state);
    fmpz_poly_init(a);
    fmpz_poly_init(b);
    fmpz_poly_init(c);
    fmpz_poly_init(d);

    fmpz_poly_randtest(a, state, 60, 1000);
    fmpz_poly_randtest(b, state, 60, 1000);

    fmpz_poly_mul(c, a, b);
    fmpz_poly_mul_classical(d, a, b);
    flint_printf("classical:       %d\n", fmpz_poly_equal(c, d));
    fmpz_poly_mul_mpn(d, a, b);
    flint_printf("mpn:             %d\n", fmpz_poly_equal(c, d));
    fmpz_poly_mul_toom_karatsuba(d, a, b);
    flint_printf("toom_karatsuba:  %d\n", fmpz_poly_equal(c, d));
    fmpz_poly_mul_KS(d, a, b);
    flint_printf("KS:              %d\n", fmpz_poly_equal(c, d));
    fmpz_poly_mul_SS(d, a, b);
    flint_printf("SS:              %d\n", fmpz_poly_equal(c, d));

    /* a window is the corresponding part of the full product */
    fmpz_poly_mulmid(d, a, b, 30, 90);
    fmpz_poly_shift_right(c, c, 30);
    fmpz_poly_truncate(c, 60);
    flint_printf("mulmid:          %d\n", fmpz_poly_equal(c, d));

    /* squaring in place */
    fmpz_poly_mul(c, a, a);
    fmpz_poly_sqr(a, a);
    flint_printf("sqr in place:    %d\n", fmpz_poly_equal(c, a));

    fmpz_poly_clear(a);
    fmpz_poly_clear(b);
    fmpz_poly_clear(c);
    fmpz_poly_clear(d);
    flint_rand_clear(state);
    return 0;
}
```
```output
classical:       1
mpn:             1
toom_karatsuba:  1
KS:              1
SS:              1
mulmid:          1
sqr in place:    1
```

## Products of many factors {#many}

To multiply many polynomials, multiply them in a balanced binary tree
(pairs, then pairs of pairs) rather than one after the other: with
fast multiplication the cost is then that of a few products of the
final size, while accumulating into one polynomial costs a quadratic
number of operations in the final size. `fmpz_poly_product_roots_fmpz_vec`
does this for linear factors. When one polynomial with large
coefficients is multiplied by many others, `fmpz_poly_mul_SS_precache_init`
keeps its transform.
