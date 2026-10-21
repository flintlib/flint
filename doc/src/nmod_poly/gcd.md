---
title: Greatest common divisors
functions: nmod_poly_gcd nmod_poly_xgcd nmod_poly_gcdinv nmod_poly_hgcd
---

Over a field $\mathbb{F}_p$, every pair of polynomials has a greatest
common divisor, unique up to a constant factor, and the non-underscore
functions return the monic one. `nmod_poly_gcd` computes $G = \gcd(A,
B)$, `nmod_poly_xgcd` also the cofactors $S, T$ with $SA + TB = G$, and
`nmod_poly_gcdinv` only the cofactor $S$ of $A$, which is what is needed
for inverses modulo $B$. The modulus must be prime; over $\Z/n\Z$ with
composite $n$ the Euclidean algorithm can meet a non-invertible leading
coefficient, which throws an exception.

::: summary algorithms/polynomials/half-gcd :::

## Algorithms {#algorithms}

The *Euclidean algorithm* replaces $(A, B)$ by $(B, A \bmod B)$ until the
remainder is zero. In the typical case each step lowers the degree by
one and costs $O(\deg B)$ operations, for $O(\deg A \deg B)$ in total.
FLINT's implementation calls `_nmod_poly_rem`, which has a dedicated
routine for quotients of length 2, and for moduli with at least half the
bits of a word (but at most $\bits - 2$ bits) and $\deg B \ge 31$ it runs
in Montgomery arithmetic (`_nmod_poly_gcd_euclidean_redc_fast`), where the
long chains of dependent multiplications are cheaper.

The *half-GCD algorithm* [@ThullYap1990] uses the fact that the quotients
of the first half of the Euclidean algorithm depend only on the top half
of the coefficients. It computes them recursively from the top
coefficients and applies them to the full polynomials as a $2 \times 2$
matrix of polynomials, so that the whole remainder sequence costs $O(M(n)
\log n)$ operations. `nmod_poly_hgcd` exposes one such step: it returns
the matrix $M$ and the two consecutive remainders $A, B$ of the sequence
that straddle half the degree of the input.

`nmod_poly_gcd` uses the half-GCD algorithm for $\mathit{lenB} \ge$
`nmod_poly_gcd_hgcd_cutoff(mod)`, between 470 and 2199 depending on the
bit size of $n$, and `nmod_poly_xgcd` for $\mathit{lenB} \ge$
`nmod_poly_xgcd_hgcd_cutoff(mod)`, between 72 and 388; the cofactors make
the Euclidean algorithm more expensive, so the extended version switches
earlier.

## Examples {#examples}

Two random polynomials of degree 2500 modulo a 50-bit prime, multiplied by
a common factor of degree 500; both algorithms find the factor, made
monic, and the cofactors of `nmod_poly_xgcd` satisfy Bézout's identity:

```c example name=nmod_poly_gcd_large
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/nmod_poly.h>

int main(void)
{
    flint_rand_t state;
    nmod_poly_t a, b, c, g, h, s, t, u;
    ulong p = n_nextprime(UWORD(1) << 49, 1);

    flint_rand_init(state);
    nmod_poly_init(a, p); nmod_poly_init(b, p); nmod_poly_init(c, p);
    nmod_poly_init(g, p); nmod_poly_init(h, p);
    nmod_poly_init(s, p); nmod_poly_init(t, p); nmod_poly_init(u, p);

    nmod_poly_rand_monic(c, state, 501);
    nmod_poly_rand(a, state, 2501);
    nmod_poly_rand(b, state, 2501);
    nmod_poly_mul(a, a, c);
    nmod_poly_mul(b, b, c);

    nmod_poly_gcd(g, a, b);               /* half-GCD at this length */
    nmod_poly_gcd_euclidean(h, a, b);
    flint_printf("deg gcd = %wd\n", nmod_poly_degree(g));
    flint_printf("hgcd == euclidean: %d\n", nmod_poly_equal(g, h));
    flint_printf("gcd == c: %d\n", nmod_poly_equal(g, c));

    nmod_poly_xgcd(h, s, t, a, b);
    nmod_poly_mul(s, s, a);
    nmod_poly_mul(t, t, b);
    nmod_poly_add(u, s, t);
    flint_printf("s a + t b == g: %d\n", nmod_poly_equal(u, h));

    nmod_poly_clear(a); nmod_poly_clear(b); nmod_poly_clear(c);
    nmod_poly_clear(g); nmod_poly_clear(h);
    nmod_poly_clear(s); nmod_poly_clear(t); nmod_poly_clear(u);
    flint_rand_clear(state);
    return 0;
}
```
```output
deg gcd = 500
hgcd == euclidean: 1
gcd == c: 1
s a + t b == g: 1
```

The common factor is found exactly because the random cofactors are
coprime, which happens with probability $1 - 1/p$.

`nmod_poly_gcdinv` returns the inverse of $A$ modulo $B$ together with
the GCD, so that one call both tests and inverts. Modulo 7, $x^2 + 1$ is
invertible modulo $x^3 + x + 1$ but not modulo $x^3 + x^2 + x + 1 = (x +
1)(x^2 + 1)$:

```c example name=nmod_poly_gcdinv
#include <flint/flint.h>
#include <flint/nmod_poly.h>

int main(void)
{
    nmod_poly_t a, b, g, s, t;

    nmod_poly_init(a, 7);
    nmod_poly_init(b, 7);
    nmod_poly_init(g, 7);
    nmod_poly_init(s, 7);
    nmod_poly_init(t, 7);

    nmod_poly_set_str(a, "3 7  1 0 1");          /* x^2 + 1 */
    nmod_poly_set_str(b, "4 7  1 1 0 1");        /* x^3 + x + 1 */
    nmod_poly_gcdinv(g, s, a, b);
    flint_printf("g = "); nmod_poly_print_pretty(g, "x");
    flint_printf(", s = "); nmod_poly_print_pretty(s, "x");
    nmod_poly_mulmod(t, s, a, b);
    flint_printf(", s a mod b = "); nmod_poly_print_pretty(t, "x");
    flint_printf("\n");

    nmod_poly_set_str(b, "4 7  1 1 1 1");        /* (x + 1)(x^2 + 1) */
    nmod_poly_gcdinv(g, s, a, b);
    flint_printf("g = "); nmod_poly_print_pretty(g, "x");
    flint_printf(", s = "); nmod_poly_print_pretty(s, "x");
    flint_printf("\n");

    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(g);
    nmod_poly_clear(s);
    nmod_poly_clear(t);
    return 0;
}
```
```output
g = 1, s = 6*x, s a mod b = 1
g = x^2+1, s = 1
```
