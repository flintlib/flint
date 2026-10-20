---
title: Proving identities
summary: Evaluating elementary expressions exactly, watching identities simplify automatically or be proved by a zero test, and verifying Machin-like formulas for π and logarithms; the example programs elementary and machin.
order: 1
---

Exact real and complex numbers ([](ca)) recognise many identities
between elementary constants. Some simplify automatically during
arithmetic, because the field in which the result is represented knows
a relation between its generators; others are recognised only by an
explicit zero test. This tutorial follows the example programs
`elementary` and `machin` (see
[](tutorials/getting-started/example-programs)).

## Elementary expressions {#elementary}

The program below evaluates three expressions. The first two simplify
during the arithmetic: $i^i$ and $e^{-\pi/2}$ are created as separate
generators, as are $\log(1 + \sqrt 2)$ and $\log(3 + 2\sqrt 2)$, but
when two of them meet in one field, the construction of the field finds
and proves the relation between them ($i^i = e^{-\pi/2}$, and
$\log(3 + 2\sqrt 2) = 2 \log(1 + \sqrt 2)$ since $3 + 2\sqrt 2 = (1 +
\sqrt 2)^2$), and the result is reduced modulo it (see
[](algorithms/exact/field-towers#ideal)). In the third expression, the field
$\mathbb{Q}(\sqrt{i}, \sqrt 2, i)$ is created with $\sqrt{i}$ as an
independent generator, and the quotient is not simplified; the printed
approximation `0e-24` only shows that the value is small. The zero test
`ca_check_is_zero` then proves that it is exactly zero, by evaluating it
as an algebraic number (see [](algorithms/exact/zero-testing)).

The second part checks Machin-like formulas $\sum_j c_j
\operatorname{atan}(1/q_j) = \pi/4$ and a formula for $\log 2$ with
hyperbolic arctangents, $\operatorname{atanh}(1/q) = \frac{1}{2}
\log((q+1)/(q-1))$. Each $\operatorname{atan}(1/q)$ is represented as
$\frac{i}{2} \log((q - i)/(q + i))$, a multiple of the logarithm of a
Gaussian rational number. When these logarithms and $\pi$ meet in one
field, the construction of the field searches for integer relations
between them with LLL and proves them exactly; the reduction modulo
these relations then gives 0 without a separate zero test.

```c example name=tutorial_exact_identities
#include <flint/calcium.h>
#include <flint/ca.h>

/* atan(1/q) */
static void
atan_inv(ca_t res, ulong q, ca_ctx_t ctx)
{
    ca_set_ui(res, q, ctx);
    ca_inv(res, res, ctx);
    ca_atan(res, res, ctx);
}

/* atanh(1/q) = log((q + 1)/(q - 1)) / 2 */
static void
atanh_inv(ca_t res, ulong q, ca_ctx_t ctx)
{
    ca_set_ui(res, q + 1, ctx);
    ca_div_ui(res, res, q - 1, ctx);
    ca_log(res, res, ctx);
    ca_div_ui(res, res, 2, ctx);
}

int main(void)
{
    ca_ctx_t ctx;
    ca_t x, y;
    slong i, j;
    /* sum c_j atan(1/q_j) = pi/4 */
    slong atan_formulas[3][4][2] = {
        {{4, 5}, {-1, 239}, {0, 0}, {0, 0}},
        {{1, 2}, {1, 5}, {1, 8}, {0, 0}},
        {{12, 49}, {32, 57}, {-5, 239}, {12, 110443}},
    };
    /* 14 atanh(1/31) + 10 atanh(1/49) + 6 atanh(1/161) = log(2) */
    slong atanh_formula[3][2] = {{14, 31}, {10, 49}, {6, 161}};

    ca_ctx_init(ctx);
    ca_init(x, ctx);
    ca_init(y, ctx);

    /* i^i - exp(-pi/2) */
    ca_i(x, ctx);
    ca_pow(x, x, x, ctx);
    ca_pi(y, ctx);
    ca_div_si(y, y, -2, ctx);
    ca_exp(y, y, ctx);
    ca_sub(x, x, y, ctx);
    flint_printf("i^i - exp(-pi/2) = "); ca_print(x, ctx); flint_printf("\n");

    /* log(1 + sqrt(2)) / log(3 + 2 sqrt(2)) */
    ca_sqrt_ui(x, 2, ctx);
    ca_add_ui(x, x, 1, ctx);
    ca_sqr(y, x, ctx);
    ca_log(x, x, ctx);
    ca_log(y, y, ctx);
    ca_div(x, x, y, ctx);
    flint_printf("log(1 + sqrt(2)) / log(3 + 2 sqrt(2)) = "); ca_print(x, ctx); flint_printf("\n");

    /* sqrt(i) sqrt(2) / (1 + i) - 1 */
    ca_i(x, ctx);
    ca_sqrt(x, x, ctx);
    ca_i(y, ctx);
    ca_add_ui(y, y, 1, ctx);
    ca_div(x, x, y, ctx);
    ca_sqrt_ui(y, 2, ctx);
    ca_mul(x, x, y, ctx);
    ca_sub_ui(x, x, 1, ctx);
    flint_printf("sqrt(i) sqrt(2) / (1 + i) - 1 = "); ca_print(x, ctx); flint_printf("\n");
    flint_printf("  is zero: "); truth_print(ca_check_is_zero(x, ctx)); flint_printf("\n");

    /* Machin-like formulas */
    for (i = 0; i < 3; i++)
    {
        ca_pi(x, ctx);
        ca_div_si(x, x, -4, ctx);
        for (j = 0; j < 4 && atan_formulas[i][j][0] != 0; j++)
        {
            atan_inv(y, atan_formulas[i][j][1], ctx);
            ca_mul_si(y, y, atan_formulas[i][j][0], ctx);
            ca_add(x, x, y, ctx);
        }
        flint_printf("atan formula %wd: ", i); ca_print(x, ctx); flint_printf("\n");
    }

    ca_set_ui(x, 2, ctx);
    ca_log(x, x, ctx);
    ca_neg(x, x, ctx);
    for (j = 0; j < 3; j++)
    {
        atanh_inv(y, atanh_formula[j][1], ctx);
        ca_mul_si(y, y, atanh_formula[j][0], ctx);
        ca_add(x, x, y, ctx);
    }
    flint_printf("atanh formula: "); ca_print(x, ctx); flint_printf("\n");

    ca_clear(x, ctx);
    ca_clear(y, ctx);
    ca_ctx_clear(ctx);
    return 0;
}
```

```output
i^i - exp(-pi/2) = 0
log(1 + sqrt(2)) / log(3 + 2 sqrt(2)) = 0.500000 {1/2}
sqrt(i) sqrt(2) / (1 + i) - 1 = 0e-24 + 0e-24*I {(-a*b*c+a*b-2)/2 where a = 0.707107 + 0.707107*I [Sqrt(1.00000*I {c})], b = 1.41421 [b^2-2=0], c = I [c^2+1=0]}
  is zero: T_TRUE
atan formula 0: 0
atan formula 1: 0
atan formula 2: 0
atanh formula: 0
```

## The example programs {#programs}

The program `elementary` evaluates about fifteen such expressions,
including $e^{\pi i} + 1$, $\log(10^{-123}) / \log 100 = -123/2$,
$2 \log(\pi i) - 4 \log \sqrt{\pi} - \pi i = 0$, an identity between
logarithms and $\pi^2$ from [@BBK2014], $\sqrt{5 + 2\sqrt 6} - \sqrt 2 -
\sqrt 3$ (proved zero by the zero test), Ramanujan's constant
$e^{\pi \sqrt{163}} - 640320^3 - 744$ (nonzero, about $-7.5 \cdot
10^{-13}$, printed as `0e-5` because `ca_print` does not refine its
approximation) and an identity between two values of the error function.
All results take a few milliseconds.

The program `machin` checks eight formulas for $\pi/4$ and seven for
$\log p$ in the same way, and then verifies the Machin-type sets
tabulated in FLINT itself, which are used for computing logarithms and
arctangents in [](mp_real): for the first $n$ primes $p_i$, relations
$D \log p_i = \sum_j C_{ij} \operatorname{atanh}(1/x_j)$, and for the
first $n$ Gaussian primes, analogous relations for their arguments.
These identities hold exactly, but `ca` proves only part of them for
larger $n$, even with the raised options `CA_OPT_LLL_PREC = 4096` and
`CA_OPT_PREC_LIMIT = 65536` that the program sets: with the default
argument 8, all rows are proved for the arctangent sets with up to six
terms and the hyperbolic sets with up to five, but only 6 of 7 and 7 of
8 rows of the larger arctangent sets, one of six rows of the six-term
hyperbolic set, and none of the rows of the seven- and eight-term hyperbolic
sets, taking about a second in total. A row that is not proved is reported as undecided, not as wrong:
the relation search does not find or cannot prove the integer relation,
with large coefficients, between the many logarithms involved, and the
zero test has no other way to prove the identity (see [](algorithms/exact/zero-testing)).
