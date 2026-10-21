---
title: Greatest common divisors
functions: nmod_mpoly_gcd nmod_mpoly_gcd_cofactors nmod_mpoly_gcd_brown nmod_mpoly_gcd_hensel nmod_mpoly_gcd_zippel nmod_mpoly_gcd_zippel2
---

Over a field $\mathbb{F}_p$, the polynomial ring $\mathbb{F}_p[x_1,
\ldots, x_m]$ is a unique factorisation domain and any two polynomials
have a GCD, unique up to a nonzero constant. `nmod_mpoly_gcd` returns the
*monic* GCD (leading coefficient 1 with respect to the ordering of the
context), and `nmod_mpoly_gcd_cofactors` also the cofactors $\bar{A} =
A/G$ and $\bar{B} = B/G$, which comes at little extra cost because all
algorithms compute them anyway. The modulus must be prime.

The functions return 1 on success. They return 0 only if the inputs
cannot be brought to exponents of at most one word per variable, even
after removing the common monomial content and dividing the exponents by
their common strides (`nmod_mpoly_deflation`); in practice this means
exponents beyond $2^{\bits}$. The outputs are then undefined.

::: summary algorithms/multivariate/modp-gcd :::

## Algorithms {#algorithms}

After the trivial cases (a zero or monomial input, cofactors that are
monomials, inputs sharing at most one variable, a variable present in
only one input), `nmod_mpoly_gcd` estimates the degrees of the GCD from
evaluations at random points and chooses, from the number of variables,
the density of the inputs and the size of $p$ compared with the degrees,
among four algorithms, trying the next one when an algorithm fails:

- `nmod_mpoly_gcd_brown`: Brown's dense algorithm [@Bro1971], which
  evaluates one variable at many points, computes the GCDs of the images
  recursively and interpolates. It is the method of choice for dense
  inputs in few variables.
- `nmod_mpoly_gcd_zippel`: Zippel's sparse interpolation [@Zip1979],
  which determines the monomials of the GCD from one image and the
  coefficients from further images by solving linear systems; efficient
  when the GCD is sparse.
- `nmod_mpoly_gcd_zippel2`: sparse interpolation with two main variables,
  using a known multiple of the leading coefficient of the GCD to fix the
  normalisation of the images [@HuMon2016].
- `nmod_mpoly_gcd_hensel`: Hensel lifting of the bivariate factorisation
  $A = G \bar{A}$, $B = G \bar{B}$ of an image to all variables (the EEZ-GCD
  approach [@Wan1980]), which needs only one good evaluation point.

When $p$ is small, random evaluation points in $\mathbb{F}_p$ often fail
or do not exist. The algorithms then work with images over extension
fields: Hensel lifting over $\mathbb{F}_{p^k}$ represented with Zech
logarithms, and Brown's and Zippel's algorithms modulo irreducible
polynomials of increasing degree in one variable, combined by the Chinese
remainder theorem. This is transparent: the result is the GCD over
$\mathbb{F}_p$ in all cases. The details, including the exact selection
rules, are described in [](algorithms/multivariate/modp-gcd).

The individual algorithms are exported mainly for testing and
benchmarking. They treat zero inputs and the trivial cases like
`nmod_mpoly_gcd`, and otherwise return 0 when the method cannot be applied
(for example when the degree bounds would overflow the dense
representation of Brown's algorithm) or fails, so a caller should be
prepared to fall back to `nmod_mpoly_gcd`.

## Example {#example}

Over $\mathbb{F}_2$, which has only two evaluation points, the GCD of two
polynomials of degree up to 15 in three variables is found by
`nmod_mpoly_gcd_cofactors` and by each of the four algorithms:

```c example name=nmod_mpoly_gcd_f2
#include <flint/flint.h>
#include <flint/nmod_mpoly.h>

int main(void)
{
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_t a, b, g, abar, bbar, t;
    const char * vars[] = {"x", "y", "z"};
    int (* algs[4])(nmod_mpoly_t, const nmod_mpoly_t, const nmod_mpoly_t,
                    const nmod_mpoly_ctx_t) = {nmod_mpoly_gcd_brown,
        nmod_mpoly_gcd_hensel, nmod_mpoly_gcd_zippel, nmod_mpoly_gcd_zippel2};
    const char * names[4] = {"brown", "hensel", "zippel", "zippel2"};
    int i;

    /* the field with two elements is far too small for random evaluation */
    nmod_mpoly_ctx_init(ctx, 3, ORD_LEX, 2);
    nmod_mpoly_init(a, ctx); nmod_mpoly_init(b, ctx); nmod_mpoly_init(g, ctx);
    nmod_mpoly_init(abar, ctx); nmod_mpoly_init(bbar, ctx); nmod_mpoly_init(t, ctx);

    nmod_mpoly_set_str_pretty(a, "(x^3*y + x*y*z + z^2 + 1)^2 * (x^4 + y*z^3 + 1)", vars, ctx);
    nmod_mpoly_set_str_pretty(b, "(x^3*y + x*y*z + z^2 + 1) * (y^5 + x*z + x)^3", vars, ctx);

    nmod_mpoly_gcd_cofactors(g, abar, bbar, a, b, ctx);
    flint_printf("gcd  = "); nmod_mpoly_print_pretty(g, vars, ctx); flint_printf("\n");
    flint_printf("abar = "); nmod_mpoly_print_pretty(abar, vars, ctx); flint_printf("\n");
    flint_printf("lengths: a %wd, b %wd, bbar %wd\n", nmod_mpoly_length(a, ctx),
                 nmod_mpoly_length(b, ctx), nmod_mpoly_length(bbar, ctx));

    for (i = 0; i < 4; i++)
    {
        int ok = algs[i](t, a, b, ctx);
        flint_printf("%-8s success %d, same gcd %d\n", names[i], ok,
                     ok && nmod_mpoly_equal(t, g, ctx));
    }

    nmod_mpoly_clear(a, ctx); nmod_mpoly_clear(b, ctx); nmod_mpoly_clear(g, ctx);
    nmod_mpoly_clear(abar, ctx); nmod_mpoly_clear(bbar, ctx); nmod_mpoly_clear(t, ctx);
    nmod_mpoly_ctx_clear(ctx);
    return 0;
}
```
```output
gcd  = x^3*y+x*y*z+z^2+1
abar = x^7*y+x^5*y*z+x^4*z^2+x^4+x^3*y^2*z^3+x^3*y+x*y^2*z^4+x*y*z+y*z^5+y*z^3+z^2+1
lengths: a 12, b 30, bbar 9
brown    success 1, same gcd 1
hensel   success 1, same gcd 1
zippel   success 1, same gcd 1
zippel2  success 1, same gcd 1
```

The cofactor $\bar{A}$ is $A/G = (x^3y + xyz + z^2 + 1)(x^4 + yz^3 + 1)$,
expanded.
