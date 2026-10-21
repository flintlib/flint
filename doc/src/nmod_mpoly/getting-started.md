---
title: Getting started with nmod_mpoly
summary: How to set up a context, build, print, inspect and free multivariate polynomials over Z/nZ, how the terms are ordered and stored, and which operations need a prime modulus.
order: 1
---

An `nmod_mpoly_t` is a polynomial in $m$ variables over $\Z/n\Z$ for a
word-size modulus $n$. Unlike `nmod_poly_t`, it does not store its
modulus: the modulus, the number of variables and the monomial ordering
belong to a context object `nmod_mpoly_ctx_t`, which is passed as the
last argument of every function. The same conventions hold
for `fmpz_mpoly`, `fmpz_mod_mpoly` and `fq_nmod_mpoly`.

## A first program {#first}

The program below works in $(\Z/7\Z)[x, y, z]$ with the degree reverse
lexicographic ordering. It reads $f$ from a string, builds $g = 3xz + y$
term by term, multiplies, inspects the product, divides it back, and
lists the terms of $f$.

```c example name=nmod_mpoly_first
#include <flint/flint.h>
#include <flint/nmod_mpoly.h>

int main(void)
{
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_t f, g, h, q;
    const char * vars[] = {"x", "y", "z"};
    ulong e[3], point[3] = {1, 2, 3};
    slong i;

    nmod_mpoly_ctx_init(ctx, 3, ORD_DEGREVLEX, 7);
    nmod_mpoly_init(f, ctx);
    nmod_mpoly_init(g, ctx);
    nmod_mpoly_init(h, ctx);
    nmod_mpoly_init(q, ctx);

    /* parse a string; integers are reduced modulo 7 */
    nmod_mpoly_set_str_pretty(f, "x^2*y + 10*y*z^3 - 1", vars, ctx);

    /* build g term by term: 3*x*z + y */
    e[0] = 1; e[1] = 0; e[2] = 1;
    nmod_mpoly_push_term_ui_ui(g, 3, e, ctx);
    e[0] = 0; e[1] = 1; e[2] = 0;
    nmod_mpoly_push_term_ui_ui(g, 1, e, ctx);
    nmod_mpoly_sort_terms(g, ctx);
    nmod_mpoly_combine_like_terms(g, ctx);

    nmod_mpoly_mul(h, f, g, ctx);

    flint_printf("f = "); nmod_mpoly_print_pretty(f, vars, ctx); flint_printf("\n");
    flint_printf("g = "); nmod_mpoly_print_pretty(g, vars, ctx); flint_printf("\n");
    flint_printf("f*g = "); nmod_mpoly_print_pretty(h, vars, ctx); flint_printf("\n");
    flint_printf("terms %wd, total degree %wd, degree in z %wd\n",
        nmod_mpoly_length(h, ctx), nmod_mpoly_total_degree_si(h, ctx),
        nmod_mpoly_degree_si(h, 2, ctx));
    flint_printf("f*g at (1,2,3) = %wu\n", nmod_mpoly_evaluate_all_ui(h, point, ctx));

    flint_printf("divides: %d, ", nmod_mpoly_divides(q, h, g, ctx));
    flint_printf("quotient == f: %d\n", nmod_mpoly_equal(q, f, ctx));

    for (i = 0; i < nmod_mpoly_length(f, ctx); i++)
    {
        nmod_mpoly_get_term_exp_ui(e, f, i, ctx);
        flint_printf("term %wd: coeff %wu, exponents (%wu, %wu, %wu)\n", i,
            nmod_mpoly_get_term_coeff_ui(f, i, ctx), e[0], e[1], e[2]);
    }

    nmod_mpoly_clear(f, ctx);
    nmod_mpoly_clear(g, ctx);
    nmod_mpoly_clear(h, ctx);
    nmod_mpoly_clear(q, ctx);
    nmod_mpoly_ctx_clear(ctx);
    return 0;
}
```
```output
f = 3*y*z^3+x^2*y+6
g = 3*x*z+y
f*g = 2*x*y*z^4+3*x^3*y*z+3*y^2*z^3+x^2*y^2+4*x*z+6*y
terms 6, total degree 6, degree in z 4
f*g at (1,2,3) = 1
divides: 1, quotient == f: 1
term 0: coeff 3, exponents (0, 1, 3)
term 1: coeff 1, exponents (2, 1, 0)
term 2: coeff 6, exponents (0, 0, 0)
```

Every polynomial is initialised with the context it will be used with,
and cleared with it before the context itself is cleared. The variable
names are only used for input and output; with `NULL` instead of `vars`,
the variables print as `x1`, `x2`, `x3`.

## Orderings and the stored form {#representation}

A polynomial is stored as an array of terms in strictly decreasing order
with respect to the monomial ordering of the context, without zero
coefficients and with coefficients in $[0, n)$. In the example, the
ordering `ORD_DEGREVLEX` compares total degrees first, so $yz^3$ (degree
4) comes before $x^2y$ (degree 3), and the leading term of $f$ is $3yz^3$.
With `ORD_LEX`, $x^2y$ would come first. Variable 0 (here $x$) is always
the most significant variable. The orderings are described in [](mpoly).

The exponent vectors are packed into fields of `bits` bits, several to a
word when they are small. The field width grows automatically when a
result needs larger exponents, so exponents are unbounded; functions that
return exponents as `ulong` or `slong` (such as
`nmod_mpoly_get_term_exp_ui` or `nmod_mpoly_degree_si`) require them to
fit, and the `fmpz` versions work in all cases.

Building a polynomial with `nmod_mpoly_set_coeff_ui_ui` keeps it sorted
but costs time linear in its length for every term. For many terms it is
better to append them with `nmod_mpoly_push_term_ui_ui`, in any order and
possibly with repeated monomials, and to call `nmod_mpoly_sort_terms` and
`nmod_mpoly_combine_like_terms` once at the end, as in the example. The
term access functions of [](nmod_mpoly#container) work with term indices
$0 \le i <$ `nmod_mpoly_length(A, ctx)`, the leading term having index 0;
those that modify exponents or coefficients can make the polynomial
non-canonical, which the same two functions repair.

## Aliasing, prime moduli and threads {#rules}

- Outputs may alias inputs (as in `nmod_mpoly_mul(f, f, g, ctx)`), except
  where a function states otherwise (for example
  `nmod_mpoly_compose_nmod_mpoly`). Several outputs of one function must
  be distinct.
- Addition, multiplication, powering, evaluation and composition work for
  every modulus $n \ge 1$. Division needs the leading coefficient of the
  divisor to be invertible modulo $n$, and throws otherwise.
  `nmod_mpoly_gcd`, `nmod_mpoly_resultant`, `nmod_mpoly_sqrt` and the
  factorisation functions of [](nmod_mpoly_factor) assume that $n$ is
  prime, which is not checked.
- Functions that return an `int` success flag fail only when intermediate
  results would need unreasonably large exponents (for example
  `nmod_mpoly_pow_fmpz` with a huge exponent and more than one term).
- `nmod_mpoly_mul`, `nmod_mpoly_divides` and `nmod_mpoly_gcd` use several
  threads for large inputs after `flint_set_num_threads` has been called
  with a value greater than 1.

The GCD functions are described in [](nmod_mpoly/gcd), and factorisation
in [](nmod_mpoly_factor/factor).
