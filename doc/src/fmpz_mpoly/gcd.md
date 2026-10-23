---
title: Greatest common divisors
functions: fmpz_mpoly_gcd fmpz_mpoly_gcd_cofactors fmpz_mpoly_gcd_brown
---

$\Z[x_0, \ldots, x_{n-1}]$ is a unique factorisation domain, so any two
polynomials have a GCD, unique up to sign. `fmpz_mpoly_gcd` returns the
GCD with positive leading coefficient (in the monomial ordering of the
context), and `fmpz_mpoly_gcd_cofactors` also returns the cofactors
$\bar A = A/G$ and $\bar B = B/G$, which several algorithms produce at no
extra cost. Both first look for structure that makes the problem easy, and
otherwise choose between modular algorithms with a cost model. The
functions `fmpz_mpoly_gcd_brown`, `fmpz_mpoly_gcd_zippel`,
`fmpz_mpoly_gcd_zippel2`, `fmpz_mpoly_gcd_hensel` and
`fmpz_mpoly_gcd_subresultant` force a particular algorithm after the same
preprocessing.

::: summary algorithms/multivariate/gcd :::

## Easy cases {#easy}

With the *monomial content* of $A$ (the GCD of its monomials) removed,
write $\operatorname{ess}(A)$ for what is left. In this order,
`fmpz_mpoly_gcd`

1. reduces inputs with exponents of more than one word per field: it
   repacks them into one word if possible, and otherwise deflates both by
   the common shifts and strides of their exponents (which does not change
   the GCD) and inflates the result;
2. handles a monomial input directly: $\gcd(c x^e, B)$ is the GCD of $c$
   with the content of $B$, times the minimum of $x^e$ and the monomial
   content of $B$;
3. checks whether $A$ and $B$ have the same number of terms and are
   proportional term by term up to monomial factors, in which case the
   cofactors are monomials;
4. divides the exponents by the GCD of their differences in each
   variable (deflation), and determines which variables occur in
   $\operatorname{ess}(A)$ and $\operatorname{ess}(B)$;
5. returns the GCD of the term contents if no variable occurs in both;
6. computes a univariate GCD with `fmpz_poly_gcd` if only one variable
   occurs in either;
7. if one of them, say $B$, is a monomial times a univariate polynomial
   $b(x_v)$ of degree at most 100, returns that monomial part times the
   GCD of $b$ with the coefficients of $A$ viewed as a polynomial in the
   other variables;
8. if a variable occurs in $\operatorname{ess}(A)$ but not in
   $\operatorname{ess}(B)$, uses $\gcd(A, B) = \gcd(\operatorname{cont}_v
   A, B)$, where $\operatorname{cont}_v A$ is the GCD of the coefficients
   of $A$ as a polynomial in that variable, recursively.

What remains are inputs in which every variable that occurs occurs in
both.

## Degree bounds and the choice of algorithm {#choice}

Next, $A$ and $B$ are reduced modulo a 63-bit prime and evaluated at
random points in all variables but one, for each variable in turn; the
degrees of the univariate GCDs give upper bounds for the degrees of $G$ in
each variable (correct with high probability) and the numbers of their
nonzero coefficients estimate the number of terms. If all bounds are zero,
the GCD is the GCD of the term contents. If the bounds equal the degrees
of $A$, the code tests whether $A$ divides $B$, and symmetrically.

Otherwise the polynomials are converted to a lex context in the variables
that occur, in an order chosen for each algorithm, and an algorithm is
selected:

- with two variables, Brown's dense modular algorithm;
- with three or more, Brown's algorithm and the Berlekamp–Massey based
  sparse algorithm (`zippel2`) are both tried, the sparse one first if its
  estimated cost is lower and the inputs are sparse (the number of
  variables times the sum of the densities of $A$ and $B$ is less than 1)
  or its estimated cost is less than 1% of Brown's;
- Zippel's sparse algorithm is the last resort.

Hensel lifting and subresultant sequences are not used by `fmpz_mpoly_gcd`;
they are available through `fmpz_mpoly_gcd_hensel` and
`fmpz_mpoly_gcd_subresultant`.

## The algorithms {#algorithms}

- **Brown** [@Bro1971] (`fmpz_mpolyl_gcd_brown`): the GCD and both
  cofactors are computed modulo word-size primes by dense evaluation and
  interpolation, variable by variable, and combined by the Chinese
  remainder theorem. The leading coefficient is fixed to the GCD $\gamma$
  of the leading coefficients of $A$ and $B$ so that the images are
  consistent. The result is proved correct without trial division: once
  the modulus $m$ exceeds twice the bound on the coefficients of $\gamma
  A$ and $\gamma B$ given by the products of the norms of $G$ and the
  cofactors, the reconstructed products agree with $\gamma A$ and $\gamma
  B$ exactly.
- **Zippel** [@Zip1979] (`fmpz_mpolyl_gcd_zippel`): after removing the
  content with respect to the main variable, the GCD modulo one prime is
  computed by recursive sparse interpolation; its support is then assumed
  for the other primes, whose images are found by solving small linear
  systems. The images are combined until they stabilise, and the result is
  verified by division. The leading coefficient problem is handled as in
  [@dKMW2005].
- **Berlekamp–Massey** [@HuMon2016] (`fmpz_mpolyl_gcd_zippel2`): with two
  main variables, the remaining variables are evaluated at successive
  powers of random values, bivariate GCDs are computed modulo a prime, and
  each coefficient of the GCD, a sparse polynomial in the remaining
  variables, is recovered from its values with the Berlekamp–Massey
  algorithm and discrete logarithms [@BenOrTiw1988]. The leading
  coefficient of the GCD in the main variables is imposed from the GCD of
  those of $A$ and $B$.
- **Hensel lifting** [@MosYun1973; @Wan1980] (`fmpz_mpolyl_gcd_hensel`):
  all variables but one are evaluated at small integers, the univariate
  GCD is computed, and the factorisation $A = G \bar A$ (or of a linear
  combination of $A$ and $B$ coprime to the cofactor) is lifted back to
  all variables.
- **Subresultants** (`fmpz_mpoly_gcd_subresultant`): the main variable is
  chosen to minimise the degrees and the number of leading and trailing
  terms, the contents are removed, and the last subresultant is made
  primitive.

All of these are described with proofs on [](algorithms/multivariate/gcd).
The modular algorithms are probabilistic in their running time (bad primes
and bad evaluation points are detected and skipped) but the results are
always correct.

## Example {#example}

$A = 6 (xy + zt - 2)^2 (x - y^3)$ and $B = -4 (xy + zt - 2) (z^2 + 3t +
1)$ have the GCD $2(xy + zt - 2)$. The integer contents 6 and 4 contribute
their GCD 2, the sign is normalised, and all algorithms agree.

```c example name=fmpz_mpoly_gcd_algorithms
#include <flint/flint.h>
#include <flint/fmpz_mpoly.h>

int main(void)
{
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t a, b, g, abar, bbar, t;
    const char * vars[] = {"x", "y", "z", "t"};
    int ok;

    fmpz_mpoly_ctx_init(ctx, 4, ORD_LEX);
    fmpz_mpoly_init(a, ctx); fmpz_mpoly_init(b, ctx); fmpz_mpoly_init(g, ctx);
    fmpz_mpoly_init(abar, ctx); fmpz_mpoly_init(bbar, ctx); fmpz_mpoly_init(t, ctx);

    /* a = 6*(x*y + z*t - 2)^2*(x - y^3), b = -4*(x*y + z*t - 2)*(z^2 + 3*t + 1) */
    fmpz_mpoly_set_str_pretty(a, "6*(x*y + z*t - 2)^2*(x - y^3)", vars, ctx);
    fmpz_mpoly_set_str_pretty(b, "-4*(x*y + z*t - 2)*(z^2 + 3*t + 1)", vars, ctx);

    ok = fmpz_mpoly_gcd_cofactors(g, abar, bbar, a, b, ctx);
    flint_printf("success: %d\n", ok);
    flint_printf("gcd  = "); fmpz_mpoly_print_pretty(g, vars, ctx); flint_printf("\n");
    flint_printf("abar = "); fmpz_mpoly_print_pretty(abar, vars, ctx); flint_printf("\n");
    flint_printf("bbar = "); fmpz_mpoly_print_pretty(bbar, vars, ctx); flint_printf("\n");

    /* the specific algorithms agree */
    {
        int (* algs[5])(fmpz_mpoly_t, const fmpz_mpoly_t, const fmpz_mpoly_t,
                        const fmpz_mpoly_ctx_t) =
            {fmpz_mpoly_gcd_brown, fmpz_mpoly_gcd_zippel, fmpz_mpoly_gcd_zippel2,
             fmpz_mpoly_gcd_hensel, fmpz_mpoly_gcd_subresultant};
        const char * names[5] = {"brown", "zippel", "zippel2", "hensel", "subresultant"};
        int i;

        for (i = 0; i < 5; i++)
        {
            ok = algs[i](t, a, b, ctx);
            flint_printf("%s: %d", names[i], ok && fmpz_mpoly_equal(t, g, ctx));
            flint_printf(i < 4 ? ", " : "\n");
        }
    }

    fmpz_mpoly_clear(a, ctx); fmpz_mpoly_clear(b, ctx); fmpz_mpoly_clear(g, ctx);
    fmpz_mpoly_clear(abar, ctx); fmpz_mpoly_clear(bbar, ctx); fmpz_mpoly_clear(t, ctx);
    fmpz_mpoly_ctx_clear(ctx);
    return 0;
}
```
```output
success: 1
gcd  = 2*x*y+2*z*t-4
abar = 3*x^2*y-3*x*y^4+3*x*z*t-6*x-3*y^3*z*t+6*y^3
bbar = -2*z^2-6*t-2
brown: 1, zippel: 1, zippel2: 1, hensel: 1, subresultant: 1
```
