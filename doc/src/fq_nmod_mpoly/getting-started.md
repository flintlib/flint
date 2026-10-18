---
title: Getting started with fq_nmod_mpoly
summary: How to set up a context over a finite field F_q, enter coefficients that involve the generator of the field, read the printed form, and factor over F_q.
order: 1
---

`fq_nmod_mpoly` works like `nmod_mpoly` (see [](nmod_mpoly/getting-started))
with coefficients in a finite field $\mathbb{F}_q = \mathbb{F}_p[a]/(f(a))$
instead of $\mathbb{F}_p$. The field is part of the polynomial context,
which is created either from an existing `fq_nmod_ctx_t`
(`fq_nmod_mpoly_ctx_init`, which copies the field, so the field context can
be cleared afterwards) or from $p$ and the degree
(`fq_nmod_mpoly_ctx_init_deg`, whose generator prints as `#`). The field
generator is available as a constant through
`fq_nmod_mpoly_set_fq_nmod_gen`; the parser `fq_nmod_mpoly_set_str_pretty`
only knows integers and the polynomial variables.

## Example {#example}

The program works over $\mathbb{F}_4 = \mathbb{F}_2[a]/(a^2 + a + 1)$.
The polynomial $x^2 + xy + y^2$ is irreducible over $\mathbb{F}_2$ but
splits over $\mathbb{F}_4$, which contains the roots $a$ and $a + 1$ of
$t^2 + t + 1$. The program then builds $(ax + y)^2$ from the generator
and takes its square root, which in characteristic 2 halves the
exponents and takes square roots of the coefficients.

```c example name=fq_nmod_mpoly_f4
#include <flint/flint.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mpoly.h>
#include <flint/fq_nmod_mpoly_factor.h>

int main(void)
{
    fq_nmod_ctx_t fctx;
    fq_nmod_mpoly_ctx_t ctx;
    fq_nmod_mpoly_t a, t;
    fq_nmod_mpoly_factor_t f;
    fq_nmod_t c;
    const char * vars[] = {"x", "y"};
    slong i;

    /* F_4 = F_2[a]/(a^2 + a + 1) */
    fq_nmod_ctx_init_ui(fctx, 2, 2, "a");
    fq_nmod_mpoly_ctx_init(ctx, 2, ORD_LEX, fctx);
    fq_nmod_ctx_clear(fctx);          /* the context keeps its own copy */

    fq_nmod_mpoly_init(a, ctx);
    fq_nmod_mpoly_init(t, ctx);
    fq_nmod_mpoly_factor_init(f, ctx);
    fq_nmod_init(c, ctx->fqctx);

    /* x^2 + x*y + y^2 is irreducible over F_2 but splits over F_4 */
    fq_nmod_mpoly_set_str_pretty(a, "x^2 + x*y + y^2", vars, ctx);
    fq_nmod_mpoly_factor(f, a, ctx);
    fq_nmod_mpoly_factor_sort(f, ctx);
    flint_printf("x^2+x*y+y^2 = ");
    fq_nmod_mpoly_factor_print_pretty(f, vars, ctx);
    flint_printf("\n");

    /* coefficients involving the generator: (a*x + y)^2 */
    fq_nmod_mpoly_set_fq_nmod_gen(t, ctx);
    fq_nmod_mpoly_gen(a, 0, ctx);
    fq_nmod_mpoly_mul(a, a, t, ctx);
    fq_nmod_mpoly_gen(t, 1, ctx);
    fq_nmod_mpoly_add(a, a, t, ctx);
    fq_nmod_mpoly_mul(a, a, a, ctx);
    flint_printf("(a*x+y)^2 = "); fq_nmod_mpoly_print_pretty(a, vars, ctx);
    flint_printf("\n");
    for (i = 0; i < fq_nmod_mpoly_length(a, ctx); i++)
    {
        fq_nmod_mpoly_get_term_coeff_fq_nmod(c, a, i, ctx);
        flint_printf("coefficient %wd: ", i); fq_nmod_print_pretty(c, ctx->fqctx);
        flint_printf("\n");
    }
    flint_printf("is square: %d\n", fq_nmod_mpoly_sqrt(t, a, ctx));
    flint_printf("root = "); fq_nmod_mpoly_print_pretty(t, vars, ctx);
    flint_printf("\n");

    fq_nmod_clear(c, ctx->fqctx);
    fq_nmod_mpoly_factor_clear(f, ctx);
    fq_nmod_mpoly_clear(a, ctx);
    fq_nmod_mpoly_clear(t, ctx);
    fq_nmod_mpoly_ctx_clear(ctx);
    return 0;
}
```
```output
x^2+x*y+y^2 = (1)
*((1)*x + (1*a)*y)^1
*((1)*x + (1*a+1)*y)^1
(a*x+y)^2 = (1*a+1)*x^2 + (1)*y^2
coefficient 0: a+1
coefficient 1: 1
is square: 1
root = (1*a)*x + (1)*y
```

The printed form of an `fq_nmod_mpoly` shows each coefficient in
parentheses as a polynomial in the generator, followed by its monomial;
`fq_nmod_mpoly_factor_print_pretty` prints the constant and then one line
per factor. Since $a^2 = a + 1$ in $\mathbb{F}_4$, $(ax + y)^2 = (a +
1)x^2 + y^2$.

## Notes {#notes}

- Coefficients are stored inline as $d$ words each ($q = p^d$), so there
  is no `fq_nmod_struct` inside a polynomial; the term access functions
  copy coefficients in and out of `fq_nmod_t` variables, which must be
  initialised with the field context `ctx->fqctx`.
- All the operations of `nmod_mpoly` that need a field (division, GCDs,
  square roots, resultants, factorisation) work over $\mathbb{F}_q$
  without further conditions. Over small fields the GCD and factorisation
  algorithms move to extensions $\mathbb{F}_{q^n}$ automatically; see
  [](algorithms/multivariate/modp-gcd) and
  [](algorithms/multivariate/modp-factoring).
- For fields small enough for Zech logarithm tables, `fq_zech_mpoly`
  offers a partial implementation with faster field arithmetic, used
  internally by the factoring code; most applications should use
  `fq_nmod_mpoly`.
