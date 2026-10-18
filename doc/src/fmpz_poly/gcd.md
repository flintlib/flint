---
title: Greatest common divisors
functions: fmpz_poly_gcd fmpz_poly_gcd_modular fmpz_poly_gcd_heuristic fmpz_poly_gcd_subresultant fmpz_poly_xgcd fmpz_poly_lcm
---

$\Z[x]$ is a unique factorisation domain, so two polynomials have a
greatest common divisor, unique up to sign. `fmpz_poly_gcd` returns the
GCD with positive leading coefficient; its content is the GCD of the
contents of the inputs, and its primitive part the GCD of their
primitive parts. `fmpz_poly_lcm` returns the least common multiple.
`fmpz_poly_xgcd` does not compute a GCD but the resultant $r$ and
cofactors with $sf + tg = r$, which is the integral form of the
extended Euclidean algorithm: over $\Z$, the GCD itself is in general
not a combination $sf + tg$ with integral $s$ and $t$.

::: summary algorithms/polynomials/gcd-integers :::

## Algorithms {#algorithms}

The Euclidean algorithm cannot be run in $\Z[x]$ directly, and running
it in $\mathbb{Q}[x]$ produces coefficients whose size grows
exponentially in the number of steps unless contents are removed at
each step. FLINT has three algorithms.

- The *subresultant algorithm* (`fmpz_poly_gcd_subresultant`) uses
  pseudo-division and divides each remainder by a known factor, so
  that the coefficients stay of the size of minors of the Sylvester
  matrix. It needs $O(n^2)$ operations on coefficients of up to
  $O(n(b + \log n))$ bits for degree $n$ and $b$-bit coefficients, and
  is used for polynomials of length less than 6.
- The *heuristic GCD* (`fmpz_poly_gcd_heuristic`) evaluates the
  primitive parts at $x = 2^k$, computes one integer GCD with GMP, and
  reads the polynomial GCD off its base-$2^k$ digits [@CGG1989]. The
  candidate is verified by checking that it divides both inputs; if it
  does not, the function returns 0. It is fast for small coefficients
  and is tried first when the coefficients of the two inputs together
  have at most 40 bits.
- The *multimodular algorithm* (`fmpz_poly_gcd_modular`) computes the
  monic GCDs of the primitive parts modulo many word-size primes with
  the fast GCD of `nmod_poly`, scales them to a common normalisation and
  combines them by the Chinese remainder theorem [@Bro1971]. Its
  running time is roughly linear in the size of the output times the
  cost of one modular GCD, and it is used in all other cases.

The result is always correct: the heuristic and multimodular
algorithms certify their output by divisibility. The multimodular
algorithm also adapts the number of primes to the size of the GCD, so
that a GCD with small coefficients, which is the common case, costs
little more than one or two modular GCDs, and a GCD equal to 1 is
usually proved by a single 20-bit prime.

The extended GCD `fmpz_poly_xgcd` (always the multimodular
`fmpz_poly_xgcd_modular`) computes the resultant first, then the
extended GCDs modulo primes that do not divide it, and stops when the
reconstructed cofactors are stable and the bound for $|sf + tg - r|$
shows that the identity holds over $\Z$.

## Examples {#examples}

GCDs keep the common content: $\gcd(6x^2 - 6, 4x + 4) = 2x + 2$. The
extended GCD of two coprime polynomials returns their resultant and
the cofactors; the LCM satisfies $fg = \pm\gcd(f, g)
\operatorname{lcm}(f, g)$.

```c example name=fmpz_poly_gcd_small
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

int main(void)
{
    fmpz_poly_t f, g, h, s, t;
    fmpz_t r;

    fmpz_poly_init(f);
    fmpz_poly_init(g);
    fmpz_poly_init(h);
    fmpz_poly_init(s);
    fmpz_poly_init(t);
    fmpz_init(r);

    fmpz_poly_set_str(f, "3  -6 0 6");        /* 6x^2 - 6 */
    fmpz_poly_set_str(g, "2  4 4");           /* 4x + 4 */
    fmpz_poly_gcd(h, f, g);
    flint_printf("gcd = "); fmpz_poly_print_pretty(h, "x"); flint_printf("\n");
    fmpz_poly_lcm(h, f, g);
    flint_printf("lcm = "); fmpz_poly_print_pretty(h, "x"); flint_printf("\n");

    fmpz_poly_set_str(f, "3  1 0 1");         /* x^2 + 1 */
    fmpz_poly_set_str(g, "3  -1 2 3");        /* 3x^2 + 2x - 1 */
    fmpz_poly_xgcd(r, s, t, f, g);
    flint_printf("r = "); fmpz_print(r);
    flint_printf(", s = "); fmpz_poly_print_pretty(s, "x");
    flint_printf(", t = "); fmpz_poly_print_pretty(t, "x");
    fmpz_poly_mul(s, s, f);
    fmpz_poly_mul(t, t, g);
    fmpz_poly_add(h, s, t);
    flint_printf(", s f + t g = "); fmpz_poly_print_pretty(h, "x");
    flint_printf("\n");

    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    fmpz_poly_clear(h);
    fmpz_poly_clear(s);
    fmpz_poly_clear(t);
    fmpz_clear(r);
    return 0;
}
```
```output
gcd = 2*x+2
lcm = 12*x^2-12
r = 20, s = 6*x+16, t = -2*x-4, s f + t g = 20
```

Two random polynomials of degree 300 with 200-bit coefficients,
multiplied by a common factor of degree 100, and, for the quadratic
subresultant algorithm, the coprime Chebyshev polynomials $T_{31}$ and
$U_{30}$ multiplied by the cyclotomic polynomial $\Phi_{35}$. All
algorithms find the common factor.

```c example name=fmpz_poly_gcd_large
#include <flint/flint.h>
#include <flint/fmpz_poly.h>

int main(void)
{
    flint_rand_t state;
    fmpz_poly_t a, b, c, g, h;

    flint_rand_init(state);
    fmpz_poly_init(a);
    fmpz_poly_init(b);
    fmpz_poly_init(c);
    fmpz_poly_init(g);
    fmpz_poly_init(h);

    fmpz_poly_randtest(a, state, 301, 200);
    fmpz_poly_randtest(b, state, 301, 200);
    fmpz_poly_randtest(c, state, 101, 50);
    fmpz_poly_set_coeff_ui(c, 100, 1);        /* degree exactly 100 */
    fmpz_poly_mul(a, a, c);
    fmpz_poly_mul(b, b, c);

    fmpz_poly_gcd(g, a, b);
    flint_printf("deg gcd = %wd\n", fmpz_poly_degree(g));
    flint_printf("gcd == c:      %d\n", fmpz_poly_equal(g, c));
    fmpz_poly_gcd_modular(h, a, b);
    flint_printf("modular:       %d\n", fmpz_poly_equal(g, h));

    fmpz_poly_chebyshev_t(a, 31);             /* coprime to U_30 */
    fmpz_poly_chebyshev_u(b, 30);
    fmpz_poly_cyclotomic(c, 35);              /* degree 24 */
    fmpz_poly_mul(a, a, c);
    fmpz_poly_mul(b, b, c);
    fmpz_poly_gcd(g, a, b);
    flint_printf("small: gcd == c: %d\n", fmpz_poly_equal(g, c));
    flint_printf("heuristic:     %d", fmpz_poly_gcd_heuristic(h, a, b));
    flint_printf(" %d\n", fmpz_poly_equal(g, h));
    fmpz_poly_gcd_subresultant(h, a, b);
    flint_printf("subresultant:  %d\n", fmpz_poly_equal(g, h));

    fmpz_poly_clear(a);
    fmpz_poly_clear(b);
    fmpz_poly_clear(c);
    fmpz_poly_clear(g);
    fmpz_poly_clear(h);
    flint_rand_clear(state);
    return 0;
}
```
```output
deg gcd = 100
gcd == c:      1
modular:       1
small: gcd == c: 1
heuristic:     1 1
subresultant:  1
```

In the first case the GCD is $c$, which is monic, because the random
cofactors are coprime and have coprime contents, which holds with high
probability.
