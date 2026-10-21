---
title: Factoring multivariate polynomials over the integers
functions: fmpz_mpoly_factor_content fmpz_mpoly_factor_squarefree fmpz_mpoly_factor fmpz_mpoly_factor_wang fmpz_mpoly_factor_zippel fmpz_mpoly_factor_zassenhaus
---

$\Z[x_0, \ldots, x_{n-1}]$ is a unique factorisation domain: every nonzero
polynomial is $\pm 1$ times a product of primes of $\Z$ and of primitive
irreducible polynomials, unique up to order and signs. `fmpz_mpoly_factor`
computes this factorisation, except that the integer content is not split
into primes: the result is the content of $A$ (with the sign of its leading
coefficient) times distinct irreducible polynomials with positive leading
coefficients, each with its multiplicity. Two coarser factorisations are
available as well, and are the first two stages of the full one.

The algorithms are described in [](algorithms/multivariate/factorisation)
and [](algorithms/multivariate/hensel-lifting); the factorisation over
$\Q$ (`fmpq_mpoly_factor`) is a thin layer over this one, and the
factorisation over finite fields is described in
[](algorithms/multivariate/modp-factoring).

## The three levels

```c example name=fmpz_mpoly_factor_levels
#include <flint/flint.h>
#include <flint/fmpz_mpoly_factor.h>

static void show(const char * name, fmpz_mpoly_factor_t f,
                 const char ** vars, fmpz_mpoly_ctx_t ctx)
{
    flint_printf("%s: ", name);
    fmpz_mpoly_factor_print_pretty(f, vars, ctx);
    flint_printf("\n");
}

int main(void)
{
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t a;
    fmpz_mpoly_factor_t f;
    const char * vars[] = {"x", "y", "z"};

    fmpz_mpoly_ctx_init(ctx, 3, ORD_LEX);
    fmpz_mpoly_init(a, ctx);
    fmpz_mpoly_factor_init(f, ctx);

    fmpz_mpoly_set_str_pretty(a,
        "-12*x^3*y^2*(y+1)*(z+1)*(x^2-y^2)^2*(x*z + y^2 + 1)", vars, ctx);

    fmpz_mpoly_factor_content(f, a, ctx);
    show("content", f, vars, ctx);
    fmpz_mpoly_factor_squarefree(f, a, ctx);
    show("squarefree", f, vars, ctx);
    fmpz_mpoly_factor(f, a, ctx);
    show("irreducible", f, vars, ctx);
    fmpz_mpoly_factor_sort(f, ctx);
    show("sorted", f, vars, ctx);
    flint_printf("matches: %d\n", fmpz_mpoly_factor_matches(a, f, ctx));

    fmpz_mpoly_factor_clear(f, ctx);
    fmpz_mpoly_clear(a, ctx);
    fmpz_mpoly_ctx_clear(ctx);
    return 0;
}
```

```output
content: -12*(x)^3*(y)^2*(z+1)^1*(y+1)^1*(x^4-2*x^2*y^2+y^4)^1*(x*z+y^2+1)^1
squarefree: -12*(x)^3*(y)^2*(z+1)^1*(y+1)^1*(x^2-y^2)^2*(x*z+y^2+1)^1
irreducible: -12*(x)^3*(y)^2*(z+1)^1*(y+1)^1*(x+y)^2*(x-y)^2*(x*z+y^2+1)^1
sorted: -12*(z+1)^1*(y+1)^1*(x*z+y^2+1)^1*(y)^2*(x-y)^2*(x+y)^2*(x)^3
matches: 1
```

`fmpz_mpoly_factor_content` removes the signed integer content and the
powers of variables, then splits the rest by repeatedly taking the content
with respect to one variable (a GCD of the coefficients of the polynomial
viewed as univariate in that variable). Here $(y+1)$ and $(z+1)$ are
contents with respect to $x$ and $(x^2 - y^2)^2$ is the content of the
remaining part with respect to $z$. The non-monomial results are
primitive with respect to each variable they involve, and pairwise
coprime, but have exponent 1 even when they are powers.

`fmpz_mpoly_factor_squarefree` then runs Yun's algorithm [@Yun1976] on each
of these pieces, with respect to a variable in which its derivative is
nonzero. In characteristic zero this gives, for each piece, the product of
its irreducible factors of each multiplicity. The bases of the result are
pairwise coprime, but one multiplicity can occur in several bases (here
$y + 1$, $z + 1$ and $xz + y^2 + 1$ all have multiplicity 1).

`fmpz_mpoly_factor` splits every squarefree base into irreducibles.
`fmpz_mpoly_factor_sort` orders the factors by exponent and then by
`fmpz_mpoly_cmp`, which makes the output canonical.

## Strategy for irreducible factors

Each squarefree primitive base $A$ is factored as follows.

1. Exponents are repacked into one word if possible (otherwise the
   function fails). The exponent vectors are then *compressed*:
   `mpoly_compression_set` looks for an integer affine change of the
   exponents that reduces the number of variables or the degrees (for
   example $x^2 y^4 + x y^2 + 1$ is a polynomial in $x y^2$), and tests the
   Newton polytope of $A$ with the criterion of Gao and Lauder
   [@GaoLau2001], which can prove $A$ irreducible from its support alone.
   After a change of variables, contents may reappear and are split off
   again. Factors found in compressed variables are mapped back with
   `fmpz_mpoly_compression_undo`.

2. In the compressed lex context, polynomials of degree 1 in the main
   variable are irreducible, and those of degree 2 are factored directly
   (a quadratic in $x_0$ factors if and only if its discriminant is a
   square, which `fmpz_mpoly_sqrt` decides). Univariate polynomials are
   factored with `fmpz_poly_factor`. Bivariate polynomials whose dense
   size $(d_0 + 1)(d_1 + 1)$ is at most 400, or that fill at least a fifth
   of that box, are factored as dense bivariate polynomials with
   `fmpz_bpoly_factor`.

3. Otherwise the other variables are evaluated at a few random small
   integers, keeping the degree in $x_0$ and squarefreeness. The degrees of
   the factors of three such univariate images are combined by
   `zassenhaus_prune`; if no subset of the factor degrees of every image
   can give a proper factor, $A$ is irreducible.

4. The leading coefficient of $A$ in $x_0$ is factored (recursively) to
   predict the leading coefficients of the factors. Then the methods are
   tried in turn, each giving up quickly when its assumptions fail:
   Wang's algorithm [@Wan1978] with evaluation at zero, if $A$ is dense
   enough (its number of terms divided by $\prod_i (d_i + 1)$ is above
   0.002) and zero was a good evaluation point in step 3; Wang's algorithm
   with random nonzero evaluation points if the density is above 0.04; the
   sparse Hensel lifting of `fmpz_mpoly_factor_irred_zippel`, which lifts
   modulo a prime with Zippel's sparse interpolation
   [@Zip1981; @Kal1985] and then lifts the coefficients $p$-adically;
   Wang's algorithm with nonzero points again; and finally
   `fmpz_mpoly_factor_irred_zassenhaus`, which lifts without leading
   coefficient information and recombines subsets of the lifted factors,
   as in the Zassenhaus algorithm for univariate polynomials.

`fmpz_mpoly_factor_wang`, `fmpz_mpoly_factor_zippel` and
`fmpz_mpoly_factor_zassenhaus` replace step 4 by a single method (Wang
with all evaluation points, sparse lifting, or recombination); steps 1–3
are the same. They give the same factorisation and exist for testing and
benchmarking. Which one is fastest depends mainly on the density of the
polynomial: dense Hensel lifting costs a polynomial in the size of the
dense representation, which is prohibitive for sparse polynomials in many
variables, while sparse lifting depends on the number of terms of the
factors [@MonTun2020].

```c example name=fmpz_mpoly_factor_methods
#include <flint/flint.h>
#include <flint/fmpz_mpoly_factor.h>

int main(void)
{
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t a;
    fmpz_mpoly_factor_t f;
    const char * vars[] = {"x", "y", "z"};
    int (* algs[4])(fmpz_mpoly_factor_t, const fmpz_mpoly_t,
                    const fmpz_mpoly_ctx_t) =
        {fmpz_mpoly_factor, fmpz_mpoly_factor_wang,
         fmpz_mpoly_factor_zippel, fmpz_mpoly_factor_zassenhaus};
    const char * names[4] = {"default", "wang", "zippel", "zassenhaus"};
    int i;

    fmpz_mpoly_ctx_init(ctx, 3, ORD_LEX);
    fmpz_mpoly_init(a, ctx);
    fmpz_mpoly_factor_init(f, ctx);

    /* the leading coefficients in x of the factors are not constant */
    fmpz_mpoly_set_str_pretty(a,
        "(2*y*z*x^2 + x + y^3 - 3*z)*(y*x^3 + z^2*x + 5)"
        "*(z*x^2 - y*z^2 + 1)^2", vars, ctx);
    flint_printf("a has %wd terms\n", fmpz_mpoly_length(a, ctx));

    for (i = 0; i < 4; i++)
    {
        int ok = algs[i](f, a, ctx);
        fmpz_mpoly_factor_sort(f, ctx);
        flint_printf("%-10s %d: ", names[i], ok);
        fmpz_mpoly_factor_print_pretty(f, vars, ctx);
        flint_printf("\n");
    }

    fmpz_mpoly_factor_clear(f, ctx);
    fmpz_mpoly_clear(a, ctx);
    fmpz_mpoly_ctx_clear(ctx);
    return 0;
}
```

```output
a has 68 terms
default    1: 1*(x^3*y+x*z^2+5)^1*(2*x^2*y*z+x+y^3-3*z)^1*(x^2*z-y*z^2+1)^2
wang       1: 1*(x^3*y+x*z^2+5)^1*(2*x^2*y*z+x+y^3-3*z)^1*(x^2*z-y*z^2+1)^2
zippel     1: 1*(x^3*y+x*z^2+5)^1*(2*x^2*y*z+x+y^3-3*z)^1*(x^2*z-y*z^2+1)^2
zassenhaus 1: 1*(x^3*y+x*z^2+5)^1*(2*x^2*y*z+x+y^3-3*z)^1*(x^2*z-y*z^2+1)^2
```

## Notes

- The functions return 0 only if a polynomial with exponents that do not
  fit one word after removing powers of variables and compression is
  met; the factorisation is then not computed. With polynomials built
  from ordinary input this does not happen.
- The algorithms are randomised (evaluation points and primes), but the
  result is deterministic: the factors are always verified, by exact
  division or by the bounds of the lifting.
- Use `fmpz_factor` on the constant to split the integer content into
  primes.
- `fmpz_mpoly_factor_squarefree` is much cheaper than `fmpz_mpoly_factor`,
  and suffices for many uses, such as computing the radical or simplifying
  exponents.
