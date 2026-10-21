---
title: Factorisation
functions: nmod_mpoly_factor nmod_mpoly_factor_squarefree nmod_mpoly_factor_separable nmod_mpoly_factor_content nmod_mpoly_factor_zassenhaus nmod_mpoly_factor_wang nmod_mpoly_factor_zippel nmod_mpoly_factor_algo
---

These functions factor a polynomial $A \in \mathbb{F}_p[x_1, \ldots,
x_m]$, $p$ prime, as
$$
A = c \prod_{i=1}^{k} f_i^{e_i},
$$
where $c \in \mathbb{F}_p$ is the leading coefficient of $A$ and the
bases $f_i$ are monic. `nmod_mpoly_factor_content` only separates the
monomial content and the contents with respect to each variable;
`nmod_mpoly_factor_squarefree` returns squarefree, pairwise coprime bases;
`nmod_mpoly_factor` returns the irreducible factorisation, in which the
bases are irreducible and pairwise distinct. The order of the factors is
not specified; `nmod_mpoly_factor_sort` sorts them by exponent and then
by `nmod_mpoly_cmp`, which makes the output canonical. The zero polynomial
gives $c = 0$ and no factors, a nonzero constant gives $c = A$ and no
factors.

All functions return 1 for success. They return 0 only in extreme cases
(exponents of more than one word after removing the content, or no usable
evaluation points or extension fields), and `f` is then undefined.

::: summary algorithms/multivariate/modp-factoring :::

## Characteristic $p$ {#char-p}

Two features of characteristic $p$ appear in the output. First, a
polynomial with zero derivative in every variable is a $p$-th power: over
$\mathbb{F}_p$, $\sum_e c_e x^{pe} = (\sum_e c_e x^e)^p$. The squarefree
decomposition therefore extracts $p$-th roots by dividing exponents, and
multiplicities divisible by $p$ come out naturally; in the example below,
$(x^3 + yz + 1)^3 = x^9 + y^3z^3 + 1$ over $\mathbb{F}_3$. Second,
factoring needs evaluation points at which the polynomial keeps its
degree and stays squarefree. When $\mathbb{F}_p$ does not have enough of
them, the polynomial is factored over an extension $\mathbb{F}_{p^k}$ and
the conjugate factors under the Frobenius automorphism are multiplied
back together. Both are invisible to the caller.

`nmod_mpoly_factor_separable` with `sep` = 1 refines the squarefree
decomposition so that each base $g$ satisfies, for each variable,
$\partial g/\partial x_i = 0$ or $\gcd(g, \partial g/\partial x_i) = 1$.
This is the input format of the irreducible factorisation.

The variants `nmod_mpoly_factor_zassenhaus`, `nmod_mpoly_factor_wang`
and `nmod_mpoly_factor_zippel` restrict the method used to lift
factorisations of bivariate images to all variables (recombination of
lifted factors, Wang's leading coefficient method with dense Hensel
lifting, or sparse Hensel lifting); `nmod_mpoly_factor_algo` takes a bit
mask of these methods. `nmod_mpoly_factor` uses all of them, choosing the
order from the density of the polynomial; see
[](algorithms/multivariate/modp-factoring#choice).

## Example {#example}

Over $\mathbb{F}_3$, the squarefree decomposition of a polynomial with a
cube factor, and its irreducible factorisation; then a factorisation over
$\mathbb{F}_2$, where $x^4 + y^4 = (x + y)^4$:

```c example name=nmod_mpoly_factor_char_p
#include <flint/flint.h>
#include <flint/nmod_mpoly.h>
#include <flint/nmod_mpoly_factor.h>

static void print_factor(const char * name, nmod_mpoly_factor_t f,
                         const char ** vars, const nmod_mpoly_ctx_t ctx)
{
    slong i;
    flint_printf("%s: %wu", name, nmod_mpoly_factor_get_constant_ui(f, ctx));
    for (i = 0; i < nmod_mpoly_factor_length(f, ctx); i++)
    {
        flint_printf(" * (");
        nmod_mpoly_print_pretty(f->poly + i, vars, ctx);
        flint_printf(")^%wd", nmod_mpoly_factor_get_exp_si(f, i, ctx));
    }
    flint_printf("\n");
}

int main(void)
{
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_t a;
    nmod_mpoly_factor_t f;
    const char * vars[] = {"x", "y", "z"};

    nmod_mpoly_ctx_init(ctx, 3, ORD_LEX, 3);
    nmod_mpoly_init(a, ctx);
    nmod_mpoly_factor_init(f, ctx);

    nmod_mpoly_set_str_pretty(a, "2*y^2*(x^3 + y*z + 1)^3*(x^2 + y^2 + z^2)*(x^2 - y^2 + z)", vars, ctx);
    flint_printf("a has %wd terms\n", nmod_mpoly_length(a, ctx));

    nmod_mpoly_factor_squarefree(f, a, ctx);
    nmod_mpoly_factor_sort(f, ctx);
    print_factor("squarefree", f, vars, ctx);

    nmod_mpoly_factor(f, a, ctx);
    nmod_mpoly_factor_sort(f, ctx);
    print_factor("irreducible", f, vars, ctx);

    /* over F_2: x^4 + y^4 = (x + y)^4, and a factorisation needing extensions */
    nmod_mpoly_ctx_set_modulus(ctx, 2);
    nmod_mpoly_set_str_pretty(a, "(x^2*y + x*z^2 + y*z + 1)*(x^3 + y^2*z + x*z + z + 1)*(x^4 + y^4)", vars, ctx);
    nmod_mpoly_factor(f, a, ctx);
    nmod_mpoly_factor_sort(f, ctx);
    print_factor("over F_2", f, vars, ctx);

    nmod_mpoly_factor_clear(f, ctx);
    nmod_mpoly_clear(a, ctx);
    nmod_mpoly_ctx_clear(ctx);
    return 0;
}
```
```output
a has 21 terms
squarefree: 2 * (x^4+x^2*z^2+x^2*z+2*y^4+2*y^2*z^2+y^2*z+z^3)^1 * (y)^2 * (x^3+y*z+1)^3
irreducible: 2 * (x^2+2*y^2+z)^1 * (x^2+y^2+z^2)^1 * (y)^2 * (x^3+y*z+1)^3
over F_2: 1 * (x^2*y+x*z^2+y*z+1)^1 * (x^3+x*z+y^2*z+z+1)^1 * (x+y)^4
```

The squarefree decomposition keeps the two factors of multiplicity 1
together as one base; the irreducible factorisation separates them. The
program switches to $\mathbb{F}_2$ with `nmod_mpoly_ctx_set_modulus`,
which is allowed because both `a` and `f` are completely overwritten
before they are used again. `nmod_mpoly_factor_get_exp_si` takes a non-`const` factorisation.
