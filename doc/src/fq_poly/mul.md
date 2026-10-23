---
title: Multiplying polynomials over finite fields
functions: fq_poly_mul fq_poly_mullow
---

`fq_poly_mul` and `fq_poly_mullow` compute the product, respectively the
truncated product, of two polynomials over $\mathbb{F}_q =
\mathbb{F}_p[y]/(f(y))$, $q = p^d$. All algorithms except the classical
one turn the problem into a single product of polynomials over $\Z$ or
over $\Z/p\Z$, so that they inherit the quasi-linear multiplication of
`fmpz_poly` and `nmod_poly` ([](algorithms/polynomials/multiplication)),
and reduce modulo $f$ only once per output coefficient.

::: summary algorithms/polynomials/kronecker :::

## Elements as polynomials {#bivariate}

An element of $\mathbb{F}_q$ is represented by a polynomial of degree
less than $d$ in the generator $y$ with coefficients in $[0, p)$, so a
polynomial of length $n$ over $\mathbb{F}_q$ is a bivariate polynomial
$A(x, y) = \sum_{i < n} \sum_{j < d} a_{ij} x^i y^j$ over $\Z$. The
product over $\mathbb{F}_q$ is obtained from the product $C = AB$ in
$\Z[x, y]$ by reducing every coefficient modulo $p$ and every coefficient
of $x^k$, a polynomial of degree at most $2d - 2$ in $y$, modulo $f$. The
coefficients of $C$ are bounded:

::: lemma {#lem-bound}
Let $A$, $B$ have lengths $n_1$, $n_2$ in $x$, degree less than $d$ in
$y$ and coefficients in $[0, p)$. Then every coefficient of $C = AB$ is
less than $\min(n_1, n_2)\, d\, (p - 1)^2 + 1$.
:::

::: proof
The coefficient of $x^k y^l$ in $C$ is $\sum a_{ij} b_{k-i, l-j}$, with
at most $\min(n_1, n_2)$ values of $i$ and at most $d$ values of $j$,
and every product is at most $(p - 1)^2$.
:::

## The algorithms {#algorithms}

**Classical** (`fq_poly_mul_classical`). Each output coefficient is a
dot product over $\mathbb{F}_q$, computed by `_fq_vec_dot` (via
`_gr_poly_mullow_classical`): the products of representatives are summed
without reduction and the sum is reduced once modulo $p$ and $f$. The
cost is $O(n_1 n_2)$ products of polynomials of degree less than $d$
plus $n_1 + n_2 - 1$ reductions.

**Kronecker substitution** (`fq_poly_mul_KS`). Every coefficient
$a_i(y)$ is packed into the integer $a_i(2^b)$ with `fq_bit_pack`, where
$b = \operatorname{bits}(d) + \operatorname{bits}(\min(n_1, n_2)) +
2 \operatorname{bits}(p)$. By [](#lem-bound) the coefficients of the
product of the packed integer polynomials, computed with
`_fmpz_poly_mul`, are the integers $c_k(2^b)$ with digits that do not
overlap, so `fq_bit_unpack` recovers $c_k(y)$ and reduces it. The
integer polynomials have length $n_1$ and $n_2$ and coefficients of about
$db$ bits; `_fmpz_poly_mul` itself uses Kronecker substitution or
multimodular FFTs, so the cost is quasi-linear in the total bit size
$n d \log p$ [@Kro1882; @GG2003, §8.4].

**Substitution $y = x^{2d-1}$** (`fq_poly_mul_univariate`). Since the
coefficients of $C$ have degree at most $2d - 2$ in $y$, the map
$y \mapsto x^{2d - 1}$ sends $A$ and $B$ to univariate polynomials whose
product determines $C$: the coefficient $c_k(y)$ occupies the block of
$2d - 1$ coefficients starting at $(2d - 1)k$. The univariate product
is computed with `_fmpz_poly_mul` for `fq` and with `_nmod_poly_mul`
(which uses `fft_small` for long products) for `fq_nmod`, and each block
is reduced modulo $f$. Compared with Kronecker substitution into
integers, no bit packing is needed and, over $\Z/p\Z$, the coefficients
stay reduced, at the price of almost doubling the length.

**Reordering** (`fq_poly_mul_reorder`). The roles of $x$ and $y$ are
exchanged: $A = \sum_{j < d} y^j A_j(x)$ with $A_j \in \Z[x]$, and the
$d^2$ products $A_j B_l$ are computed with fast integer polynomial
multiplication and combined classically in the $y$ direction before the
reduction modulo $f$. This pays off for long polynomials over fields of
small degree $d$.

**Zech logarithms.** In `fq_zech_poly` the coefficients are discrete
logarithms, so the classical algorithm uses table lookups only. The
substitution algorithm first converts the coefficients to their
polynomial representation with a table (`eval_table`), multiplies with
`_nmod_poly_mul`, reduces the blocks modulo $f$ with word arithmetic
(accumulating lazily when $(d - 1)(p - 1)^2 + p$ fits in a word), and
converts the results back with $d$ Zech table lookups per coefficient.

## Algorithm selection {#selection}

| Module | `mul` | `mullow` | `sqr` |
|---|---|---|---|
| `fq_poly` | classical if both lengths are below 6, reordering if $d < 4$, Kronecker substitution otherwise | classical if $n < 6$ or both lengths are below 6, Kronecker substitution otherwise | as `mul` |
| `fq_nmod_poly` | classical if a length is 1, if $d = 2$ and both lengths are at most 2, or if $\operatorname{bits}(p) \cdot d \cdot \max(n_1, n_2) \le 8$; substitution $y = x^{2d-1}$ otherwise | classical if $n < 6$ or both lengths are below 6, Kronecker substitution otherwise | classical below length 6, Kronecker substitution otherwise |
| `fq_zech_poly` | classical if a length is below $2d + 2$; otherwise substitution $y = x^{2d-1}$ if `_fq_zech_poly_mulmid_want_univariate` predicts it to be faster, classical otherwise | as `mul`, with the work estimate for the truncated product | as `mul`, with `_fq_zech_poly_sqr_want_univariate` |

The threshold of `fq_zech_poly` (`_fq_zech_poly_mul_univariate_threshold`)
is $9d - 4$ terms for $q < 2^{17}$ and decreases to $2d + 2$ for larger
fields, where the tables no longer fit in the cache and lookups become
expensive. `fq_poly_mulhigh` always uses the classical algorithm, while
`_fq_poly_mulhigh` computes a full product by Kronecker substitution for
inputs of length at least 6.

## Example {#example}

All algorithms give the same product. The program multiplies random
polynomials of lengths 40 and 30 over $\mathbb{F}_{101^3}$ with each
algorithm:

```c example name=fq_poly_mul_algorithms
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fq.h>
#include <flint/fq_poly.h>

int main(void)
{
    fmpz_t p;
    fq_ctx_t ctx;
    fq_poly_t a, b, c, d;
    flint_rand_t state;

    fmpz_init_set_ui(p, 101);
    fq_ctx_init(ctx, p, 3, "y");               /* GF(101^3) */
    flint_rand_init(state);
    fq_poly_init(a, ctx);
    fq_poly_init(b, ctx);
    fq_poly_init(c, ctx);
    fq_poly_init(d, ctx);

    fq_poly_randtest(a, state, 40, ctx);
    fq_poly_randtest(b, state, 30, ctx);

    fq_poly_mul(c, a, b, ctx);
    fq_poly_mul_classical(d, a, b, ctx);
    flint_printf("classical:  %d\n", fq_poly_equal(c, d, ctx));
    fq_poly_mul_KS(d, a, b, ctx);
    flint_printf("KS:         %d\n", fq_poly_equal(c, d, ctx));
    fq_poly_mul_univariate(d, a, b, ctx);
    flint_printf("univariate: %d\n", fq_poly_equal(c, d, ctx));
    fq_poly_mul_reorder(d, a, b, ctx);
    flint_printf("reorder:    %d\n", fq_poly_equal(c, d, ctx));

    fq_poly_mullow(d, a, b, 20, ctx);
    fq_poly_truncate(c, 20, ctx);
    flint_printf("mullow:     %d\n", fq_poly_equal(c, d, ctx));

    fq_poly_clear(a, ctx);
    fq_poly_clear(b, ctx);
    fq_poly_clear(c, ctx);
    fq_poly_clear(d, ctx);
    flint_rand_clear(state);
    fq_ctx_clear(ctx);
    fmpz_clear(p);
    return 0;
}
```
```output
classical:  1
KS:         1
univariate: 1
reorder:    1
mullow:     1
```
