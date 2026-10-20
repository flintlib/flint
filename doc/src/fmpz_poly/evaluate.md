---
title: Evaluation
functions: fmpz_poly_evaluate_fmpz fmpz_poly_evaluate_divconquer_fmpz fmpz_poly_evaluate_fmpq fmpz_poly_evaluate_divconquer_fmpq fmpz_poly_evaluate_mod fmpz_poly_evaluate_fmpz_vec
---

`fmpz_poly_evaluate_fmpz` computes $f(a)$ for an integer $a$,
`fmpz_poly_evaluate_fmpq` for a rational $a$ (in lowest terms), and
`fmpz_poly_evaluate_mod` the value modulo a word-size integer. Values
at many points are computed one at a time by
`fmpz_poly_evaluate_fmpz_vec`. For numerical evaluation with error
bounds at real or complex points, use `arb_fmpz_poly_evaluate_arb` and
`arb_fmpz_poly_evaluate_acb`.

## Horner's rule and divide and conquer {#algorithms}

Horner's rule $f(a) = (\cdots(f_n a + f_{n-1}) a + \cdots) a + f_0$
uses $n$ multiplications, but the intermediate results grow: for $a$
of $k$ bits, the $i$-th product multiplies a number of about $ik$ bits
by one of $k$ bits, so the cost is $O(n^2 k^2)$ with classical integer
multiplication and still quadratic in $n$ with fast multiplication,
since one operand is always small. The divide-and-conquer algorithm
writes $f = f_0 + x^m f_1$ with $m = n/2$, evaluates the halves
recursively and returns $f_0(a) + a^m f_1(a)$, with the powers
$a^{2^j}$ computed once. Its products are balanced, so with fast
integer multiplication the total cost is $O(M(nk) \log n)$ for the
cost $M$ of multiplying integers of $nk$ bits. `fmpz_poly_evaluate_fmpz`
switches from Horner's rule at length 50.

At a rational point $a = p/q$, both algorithms compute the integer
$q^n f(p/q) = \sum_i f_i p^i q^{n-i}$ and divide out the common factor
at the end. `fmpz_poly_evaluate_fmpq` uses Horner's rule for length
less than 40 or when the denominator has more than $0.003 n^2$ bits,
and divide and conquer otherwise.

Evaluation modulo $n$ reduces each coefficient once and runs Horner's
rule in word arithmetic. To evaluate a polynomial modulo $n$ at many
points, convert it with `fmpz_poly_get_nmod_poly` and use
`nmod_poly_evaluate_nmod_vec`, which is quasi-linear.

## Examples {#examples}

```c example name=fmpz_poly_evaluate_points
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>

int main(void)
{
    fmpz_poly_t f;
    fmpz_t a, y;
    fmpq_t b, z;

    fmpz_poly_init(f);
    fmpz_init(a);
    fmpz_init(y);
    fmpq_init(b);
    fmpq_init(z);

    fmpz_poly_set_str(f, "4  -1 -1 0 1");     /* x^3 - x - 1 */

    fmpz_set_si(a, -3);
    fmpz_poly_evaluate_fmpz(y, f, a);
    flint_printf("f(-3) = "); fmpz_print(y); flint_printf("\n");

    fmpq_set_si(b, 4, 3);
    fmpz_poly_evaluate_fmpq(z, f, b);
    flint_printf("f(4/3) = "); fmpq_print(z); flint_printf("\n");

    flint_printf("f(10^6) mod 1000003 = %wu\n",
                 fmpz_poly_evaluate_mod(f, 1000000, 1000003));

    fmpz_poly_clear(f);
    fmpz_clear(a);
    fmpz_clear(y);
    fmpq_clear(b);
    fmpq_clear(z);
    return 0;
}
```
```output
f(-3) = -25
f(4/3) = 1/27
f(10^6) mod 1000003 = 999978
```

The divide-and-conquer algorithm and Horner's rule give the same
value. The value of a polynomial of degree 999 with 10-bit
coefficients at $3^{100}$ has about 160000 bits:

```c example name=fmpz_poly_evaluate_large
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

int main(void)
{
    flint_rand_t state;
    fmpz_poly_t f;
    fmpz_t a, y1, y2;

    flint_rand_init(state);
    fmpz_poly_init(f);
    fmpz_init(a);
    fmpz_init(y1);
    fmpz_init(y2);

    fmpz_poly_randtest(f, state, 1000, 10);
    fmpz_poly_set_coeff_ui(f, 999, 1);
    fmpz_ui_pow_ui(a, 3, 100);

    fmpz_poly_evaluate_horner_fmpz(y1, f, a);
    fmpz_poly_evaluate_divconquer_fmpz(y2, f, a);
    flint_printf("equal: %d, bits: %wd\n", fmpz_equal(y1, y2), fmpz_bits(y1));

    fmpz_poly_clear(f);
    fmpz_clear(a);
    fmpz_clear(y1);
    fmpz_clear(y2);
    flint_rand_clear(state);
    return 0;
}
```
```output
equal: 1, bits: 158338
```
