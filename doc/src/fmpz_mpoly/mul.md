---
title: Multiplication
functions: fmpz_mpoly_mul fmpz_mpoly_mul_johnson fmpz_mpoly_mul_array fmpz_mpoly_mul_dense fmpz_mpoly_mul_monomial
---

`fmpz_mpoly_mul` is the function to use for products in $\Z[x_0, \ldots,
x_{n-1}]$. Depending on the inputs it uses Johnson's heap method, which
is efficient for sparse polynomials and needs memory only for the inputs
and the output; a dense array method, which accumulates the products in a
dense array one slice at a time and is several times faster when the
product is reasonably dense; or a Kronecker substitution into a
univariate product, which is fastest for dense inputs of moderate degree.
Large products are computed with several threads. The other functions
documented below run one of these algorithms directly, which is useful for
benchmarking and for code that knows its inputs.

::: summary algorithms/multivariate/multiplication :::

## Algorithm selection {#selection}

Let $B$ and $C$ have $\ell_B$ and $\ell_C$ terms. `fmpz_mpoly_mul`
proceeds as follows.

1. If $B$ or $C$ is zero, the product is zero; if one of them has a
   single term, `fmpz_mpoly_mul_monomial` multiplies each term of the
   other one by it.
2. With one variable and single word exponents, the product is computed
   as a dense univariate product with `_fmpz_poly_mul`, provided its degree
   is at most $\ell_B \ell_C$ (and at most $\ell_B \ell_C/4$ if the
   coefficients of $B$ and $C$ together have more than `FLINT_BITS` bits).
3. If $\min(\ell_B, \ell_C) < 20$ or $\max(\ell_B, \ell_C) < 50$, the
   single-threaded heap method is used.
4. If an exponent field needs more than one word, the heap method is used,
   with threads for large inputs.
5. Otherwise, with $D = \prod_i (\deg_{x_i} B + \deg_{x_i} C + 1)$, the
   size of the dense box containing the product:
   - the *array method* is considered if $2 \le n < 8$ and the packed
     monomials of both inputs fit one word: for `ORD_LEX` if $D \le 5
     \cdot 10^7$ and $D/(\ell_B \ell_C) < 10$, and for the degree orderings
     if the number of monomials of degree at most the total degree
     $\deg B + \deg C$ of the product, estimated as
     $(\deg B + \deg C)^n/n!$, is at most $5 \cdot 10^6$
     and less than $10 \ell_B \ell_C$;
   - the *dense method* (`fmpz_mpoly_mul_dense`) is tried first if $D <
     2^{\min(\bits/3 + 16, \bits - 3)}$ and $D < \ell_B \ell_C/128$ if
     the array method is possible, or $D < \ell_B \ell_C/32$ if not;
   - otherwise the array method is used if it was found possible, and the
     heap method if not or if the array method declines.

The ratios encode the assumption that the dense method costs about $D$
operations, the heap method about $\ell_B \ell_C$ heap operations, and
the array method about four times less than the heap method. Threads are
requested for the array and heap methods as `flint_request_threads(...,
min(lB, lC)/512)`, that is, only when both inputs have more than 512
terms; whether threads are available is set by `flint_set_num_threads`.

The heap method needs memory proportional to $\ell_B + \ell_C$ plus the
output. The array method needs a dense array of at most 300000
coefficients for one slice (the main variable in lex, or the total degree)
of the product, times one to three words per coefficient. The dense method
needs memory proportional to $D$.

## The algorithms {#algorithms}

*Johnson's heap method* [@Joh1974] generates the products $b_i c_j$ in
decreasing order of their monomials by keeping, for each term $b_i$ of
the shorter input, the next product $b_i c_j$ not yet used in a heap.
Equal monomials are merged in the heap (*chaining*), coefficients of
products of word-size coefficients are accumulated in three words without
any `fmpz` arithmetic, and the exponents are compared as packed words, as
proposed by Monagan and Pearce [@MonPea2007]. The parallel version
[@MonPea2009] splits the output into ranges of monomials of decreasing
size, computed by worker threads that each run a heap over the products
in their range, and concatenates the results.

The *array method* views the product as a polynomial in the main variable
(lex) or graded by total degree (degree orderings). For each power of the
main variable, respectively each total degree, of the product it adds all
contributing products of slices of $B$ and $C$ into a dense array indexed
by the remaining exponents, then reads off the nonzero entries in order.
Coefficients are accumulated in one, two or three words when a bound on
their size allows it, and as `fmpz` otherwise. Threads compute different
slices.

The *dense method* is Kronecker substitution: with $r_i = \deg_{x_i} B +
\deg_{x_i} C + 1$, the substitution $x_0 \mapsto x$, $x_1 \mapsto
x^{r_0}$, $x_2 \mapsto x^{r_0 r_1}$, ... maps $B$ and $C$ to univariate
polynomials whose product determines $BC$ uniquely; it is computed with
`fmpz_poly_mul`.

## Examples {#examples}

The first product is a scaled-down version of Fateman's benchmark: $f =
(1 + x + y + z + t)^{10}$, which has 1001 terms, times $f + 1$. All
methods apply and agree. The second product is sparse with large degrees,
where only the heap method applies.

```c example name=fmpz_mpoly_mul_methods
#include <flint/flint.h>
#include <flint/fmpz_mpoly.h>

int main(void)
{
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t f, g, p1, p2, p3, p4;
    const char * vars[] = {"x", "y", "z", "t"};
    int ok3, ok4;

    fmpz_mpoly_ctx_init(ctx, 4, ORD_DEGREVLEX);
    fmpz_mpoly_init(f, ctx); fmpz_mpoly_init(g, ctx);
    fmpz_mpoly_init(p1, ctx); fmpz_mpoly_init(p2, ctx);
    fmpz_mpoly_init(p3, ctx); fmpz_mpoly_init(p4, ctx);

    /* Fateman's benchmark, scaled down: f = (1 + x + y + z + t)^10 */
    fmpz_mpoly_set_str_pretty(f, "(1 + x + y + z + t)^10", vars, ctx);
    fmpz_mpoly_add_ui(g, f, 1, ctx);

    fmpz_mpoly_mul(p1, f, g, ctx);             /* automatic choice */
    fmpz_mpoly_mul_johnson(p2, f, g, ctx);     /* heap */
    ok3 = fmpz_mpoly_mul_array(p3, f, g, ctx); /* dense array, may decline */
    ok4 = fmpz_mpoly_mul_dense(p4, f, g, ctx); /* Kronecker substitution */

    flint_printf("f has %wd terms, f*g has %wd terms\n",
                 fmpz_mpoly_length(f, ctx), fmpz_mpoly_length(p1, ctx));
    flint_printf("johnson == mul: %d\n", fmpz_mpoly_equal(p2, p1, ctx));
    flint_printf("array: %d, equal: %d\n", ok3, ok3 && fmpz_mpoly_equal(p3, p1, ctx));
    flint_printf("dense: %d, equal: %d\n", ok4, ok4 && fmpz_mpoly_equal(p4, p1, ctx));

    /* a sparse product with large degrees: only the heap method applies */
    fmpz_mpoly_set_str_pretty(f, "x^1000000 + y^1000000 + z^1000000 + t^1000000 + 1", vars, ctx);
    fmpz_mpoly_set_str_pretty(g, "x^999999*y + y^999999*z + z^999999*t + 1", vars, ctx);
    fmpz_mpoly_mul(p1, f, g, ctx);
    flint_printf("sparse: %wd terms, array: %d, dense: %d\n", fmpz_mpoly_length(p1, ctx),
                 fmpz_mpoly_mul_array(p3, f, g, ctx), fmpz_mpoly_mul_dense(p4, f, g, ctx));

    fmpz_mpoly_clear(f, ctx); fmpz_mpoly_clear(g, ctx);
    fmpz_mpoly_clear(p1, ctx); fmpz_mpoly_clear(p2, ctx);
    fmpz_mpoly_clear(p3, ctx); fmpz_mpoly_clear(p4, ctx);
    fmpz_mpoly_ctx_clear(ctx);
    return 0;
}
```
```output
f has 1001 terms, f*g has 10626 terms
johnson == mul: 1
array: 1, equal: 1
dense: 1, equal: 1
sparse: 20 terms, array: 0, dense: 0
```

When called directly, `fmpz_mpoly_mul_heap_threaded` and
`fmpz_mpoly_mul_array_threaded` request threads as
`flint_request_threads(..., min(lB, lC)/16)`.
