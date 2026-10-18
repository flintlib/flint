---
title: Factoring polynomials over a prime field
functions: nmod_poly_factor nmod_poly_factor_with_cantor_zassenhaus nmod_poly_factor_cantor_zassenhaus
---

`nmod_poly_factor` computes the complete factorisation
$$
f = c \prod_{i=1}^{s} f_i^{\,e_i}
$$
of a polynomial $f \in \mathbb{F}_p[x]$ into distinct monic irreducible
polynomials $f_i$ with multiplicities $e_i \ge 1$: the factors and
exponents are stored in an `nmod_poly_factor_t`, and the leading
coefficient $c$ is the return value. The modulus $p$ must be prime. The
`_with_` variants force one algorithm and are otherwise identical; the
three `void` functions append to an existing factorisation instead of
overwriting it and discard $c$, which is convenient when factoring the
parts of a product one by one.

::: summary algorithms/polynomials/factoring-finite-fields :::

## Algorithms {#algorithms}

All functions are wrappers of the generic implementation in `gr_poly`
over the `nmod` context, which proceeds in four steps:

1. *Deflation.* If $f(x) = g(x^m)$ with $m \ge 2$, $g$ is factored first,
   and then each $g_j(x^m)$ separately.
2. *Squarefree factorisation*, by Yun-type gcd iterations and $p$-th roots
   ([](algorithms/polynomials/squarefree)), which costs a few gcds.
3. *Distinct-degree factorisation* of each squarefree part, which
   separates the irreducible factors by degree.
4. *Equal-degree factorisation*, a randomised splitting of each product of
   irreducible factors of the same degree.

The algorithms differ in steps 3 and 4. The Cantor–Zassenhaus algorithm
finds the factors of degree $d$ as $\gcd(g, x^{p^d} - x)$ for $d = 1, 2,
\ldots$, computing each $x^{p^d} \bmod g$ by an exponentiation by $p$; its
cost is about $\tfrac{n}{2} \cdot \tfrac{3}{2} \log_2 p$ multiplications
modulo $g$ for an irreducible $g$ of degree $n$, so it suits small $p$ and
small degrees. The Kaltofen–Shoup algorithm replaces the exponentiations
by modular compositions in a baby-step giant-step scheme, about
$\sqrt{2n}$ compositions and $n/2$ multiplications modulo $g$, independent
of $p$. Berlekamp's algorithm computes the kernel of an $n \times n$
matrix over $\mathbb{F}_p$ instead; at $O(n^3)$ (or $O(n^\omega)$) it is
slower for large $n$, and is mainly useful for checking the others.

`nmod_poly_factor` chooses Cantor–Zassenhaus when the degree $n$ of the
(deflated) input satisfies $n < 125 / (\log_2 p)^{1.3}$ and
Kaltofen–Shoup otherwise:

| $p$ | 2 | 3 | 7 | 31 | $\approx 2^8$ | $\approx 2^{16}$ | $\ge 2^{40}$ |
|---|---|---|---|---|---|---|---|
| Cantor–Zassenhaus for $n <$ | 125 | 69 | 33 | 16 | 9 | 4 | 2 |

For large inputs the work is dominated by the distinct-degree
factorisation, which uses several threads when they are available
(`flint_set_num_threads`); see
[](algorithms/polynomials/factoring-finite-fields#choice).

The factors are returned in no particular order: grouped by squarefree
part, and within a part in the order in which the distinct-degree and
equal-degree steps find them. The random choices use a random state with
a fixed seed, so repeated calls give the same output.

## Examples {#examples}

A product of known factors modulo 7, one of them with multiplicity $7 =
p$, is recovered together with its leading coefficient. The squarefree
factorisation alone separates the multiplicities:

```c example name=nmod_poly_factor_basic
#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>

int main(void)
{
    nmod_poly_t f, g;
    nmod_poly_factor_t fac;
    ulong c;

    nmod_poly_init(f, 7);
    nmod_poly_init(g, 7);
    nmod_poly_factor_init(fac);

    /* f = 3 x (x^2 + 1)^2 (x + 1)^7 (x^3 + x + 1) */
    nmod_poly_set_str(f, "3 7  1 0 1");
    nmod_poly_pow(f, f, 2);
    nmod_poly_set_str(g, "2 7  1 1");
    nmod_poly_pow(g, g, 7);
    nmod_poly_mul(f, f, g);
    nmod_poly_set_str(g, "4 7  1 1 0 1");
    nmod_poly_mul(f, f, g);
    nmod_poly_shift_left(f, f, 1);
    nmod_poly_scalar_mul_nmod(f, f, 3);

    c = nmod_poly_factor(fac, f);
    flint_printf("leading coefficient %wu\n", c);
    nmod_poly_factor_print_pretty(fac, "x");

    nmod_poly_factor_squarefree(fac, f);
    flint_printf("squarefree:\n");
    nmod_poly_factor_print_pretty(fac, "x");

    nmod_poly_clear(f);
    nmod_poly_clear(g);
    nmod_poly_factor_clear(fac);
    return 0;
}
```
```output
leading coefficient 3
x ^ 1
x^3+x+1 ^ 1
x^2+1 ^ 2
x+1 ^ 7
squarefree:
x^4+x^2+x ^ 1
x^2+1 ^ 2
x+1 ^ 7
```

The factor $(x + 1)^7$ has a multiplicity divisible by $p = 7$; the
squarefree factorisation finds it by taking a $p$-th root.

By [](algorithms/polynomials/factoring-finite-fields#thm-xqd),
$x^{3^6} - x$ is the product of all monic irreducible polynomials over
$\mathbb{F}_3$ whose degree divides 6. Factoring it and counting the
factors by degree gives the numbers $I_d = \frac1d \sum_{e \mid d}
\mu(d/e)\, 3^e$ of monic irreducible polynomials of degree $d$; all three
algorithms agree:

```c example name=nmod_poly_factor_count
#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>

static void count(const nmod_poly_factor_t fac)
{
    slong i, num[7] = {0};
    for (i = 0; i < fac->num; i++)
        num[nmod_poly_degree(fac->p + i)]++;
    flint_printf("%wd %wd %wd %wd %wd %wd\n",
        num[1], num[2], num[3], num[4], num[5], num[6]);
}

int main(void)
{
    nmod_poly_t f;
    nmod_poly_factor_t fac;

    nmod_poly_init(f, 3);
    nmod_poly_factor_init(fac);

    nmod_poly_set_coeff_ui(f, 729, 1);
    nmod_poly_set_coeff_ui(f, 1, 2);       /* x^729 - x */

    nmod_poly_factor(fac, f);
    count(fac);
    nmod_poly_factor_with_cantor_zassenhaus(fac, f);
    count(fac);
    nmod_poly_factor_with_berlekamp(fac, f);
    count(fac);

    nmod_poly_clear(f);
    nmod_poly_factor_clear(fac);
    return 0;
}
```
```output
3 3 8 0 0 116
3 3 8 0 0 116
3 3 8 0 0 116
```
