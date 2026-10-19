---
title: Working with power series rings
summary: A short program that computes in Q[[x]] with error terms and in Q[[x]]/(x^n), and the difference between the two rings for zero tests.
order: 1
---

The two power series rings of `gr_series`, $R[[x]]$ with error terms
and the truncated ring $R[[x]]/(x^n)$, are generic rings, so their
elements are created with `gr_heap_init` (or `GR_TMP_INIT`) and
manipulated with the functions of `gr.h` and `gr_special.h`. The
program below works in $\mathbb{Q}[[x]]$ with default precision 6 and
in $\mathbb{Q}[[x]]/(x^6)$.

```c example name=gr_series_basics
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_poly.h>
#include <flint/gr_series.h>
#include <flint/gr_special.h>

int main(void)
{
    gr_ctx_t QQ, S, M;
    gr_ptr x, f, g;
    truth_t t;

    gr_ctx_init_fmpq(QQ);

    /* power series Q[[x]] with default precision 6 */
    gr_series_ctx_init(S, QQ, 6);
    x = gr_heap_init(S);
    f = gr_heap_init(S);
    g = gr_heap_init(S);

    GR_MUST_SUCCEED(gr_gen(x, S));
    GR_MUST_SUCCEED(gr_exp(f, x, S));
    flint_printf("exp(x)       = "); gr_println(f, S);
    GR_MUST_SUCCEED(gr_sub_ui(g, x, 1, S));
    GR_MUST_SUCCEED(gr_inv(g, g, S));
    flint_printf("1/(x-1)      = "); gr_println(g, S);
    GR_MUST_SUCCEED(gr_mul(g, x, x, S));
    flint_printf("x*x          = "); gr_println(g, S);
    GR_MUST_SUCCEED(gr_sub(g, f, f, S));
    flint_printf("f - f        = "); gr_println(g, S);
    t = gr_is_zero(g, S);
    flint_printf("f - f == 0 ? %s\n", t == T_TRUE ? "true" : t == T_FALSE ? "false" : "unknown");

    gr_heap_clear(x, S);
    gr_heap_clear(f, S);
    gr_heap_clear(g, S);
    gr_ctx_clear(S);

    /* truncated power series Q[[x]]/(x^6) */
    gr_series_mod_ctx_init(M, QQ, 6);
    x = gr_heap_init(M);
    f = gr_heap_init(M);
    GR_MUST_SUCCEED(gr_gen(x, M));
    GR_MUST_SUCCEED(gr_exp(f, x, M));
    GR_MUST_SUCCEED(gr_sub(f, f, f, M));
    flint_printf("f - f mod x^6 = "); gr_println(f, M);
    t = gr_is_zero(f, M);
    flint_printf("f - f == 0 ? %s\n", t == T_TRUE ? "true" : t == T_FALSE ? "false" : "unknown");
    gr_heap_clear(x, M);
    gr_heap_clear(f, M);
    gr_ctx_clear(M);

    gr_ctx_clear(QQ);
    return 0;
}
```
```output
exp(x)       = 1 + x + (1/2)*x^2 + (1/6)*x^3 + (1/24)*x^4 + (1/120)*x^5 + O(x^6)
1/(x-1)      = -1 - x - x^2 - x^3 - x^4 - x^5 + O(x^6)
x*x          = x^2
f - f        = 0 + O(x^6)
f - f == 0 ? unknown
f - f mod x^6 = 0 (mod x^6)
f - f == 0 ? true
```

The generator $x$ and the product $x \cdot x$ are exact series: $x^2$
has no error term because it is shorter than the default precision. The
exponential and the inverse are infinite series and carry the error
term $O(x^6)$. Subtracting $\exp(x)$ from itself gives $0 + O(x^6)$,
which is not known to be zero: the two operands are the same element,
but in general the tails hidden by two error terms $O(x^6)$ need not
cancel. In the truncated ring, where $x^6 = 0$ by definition, the same
computation gives an exact zero.

Over an inexact base ring such as `arb`, the coefficients are balls and
the error term covers the tail, so a `gr_series` element over `arb`
encloses a set of power series in both senses; the special functions
of `gr_series` (gamma functions, zeta functions, hypergeometric
functions and others) are available over `arb` and `acb`.
