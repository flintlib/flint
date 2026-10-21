---
title: Roots modulo a prime and modulo composite numbers
functions: nmod_poly_roots nmod_poly_roots_factored
---

`nmod_poly_roots` finds all roots of a nonzero polynomial $f$ in the
prime field $\mathbb{F}_p$ and returns them as the linear factors $x - a$
of an `nmod_poly_factor_t`, with the multiplicity of each root as its
exponent if requested. `nmod_poly_roots_factored` does the same over
$\Z/n\Z$ for a composite modulus $n$ whose prime factorisation is known.
Over $\Z/n\Z$ the number of roots is not bounded by the degree, and the
multiplicity of $a$ is defined as the largest $e$ with $(x - a)^e \mid f$
in $(\Z/n\Z)[x]$.

To check whether a polynomial splits into distinct linear factors with
nonzero roots, and find them if it does, the specialised
`nmod_poly_find_distinct_nonzero_roots` in `nmod_poly` is faster, since it
skips the squarefree factorisation (see [](nmod_poly/roots)).

::: summary algorithms/polynomials/root-finding :::

## Algorithms {#algorithms}

Modulo a prime, `nmod_poly_roots` computes the squarefree factorisation
of $f$, then for each squarefree part $g$ the two gcds $\gcd(g,
x^{(p-1)/2} - 1)$ and $\gcd(g, x^{(p-1)/2} + 1)$, whose product is the
product of the distinct nonzero linear factors of $g$, and splits them
with random shifts $\gcd(u, (x + c)^{(p-1)/2} - 1)$. The cost is
dominated by the exponentiations, $O(M(n) \log p \log n)$ in expectation
for degree $n$; the irreducible factors of degree at least 2 are never
separated from each other. The root 0 is split off directly.

Modulo $n = \prod_j p_j^{k_j}$, `nmod_poly_roots_factored` finds the roots
modulo each $p_j$, lifts every one of them to all roots modulo
$p_j^{k_j}$ that reduce to it, using
[](algorithms/polynomials/root-finding#lem-hensel) with doubling
precision, and combines the roots for different primes with the Chinese
remainder theorem. A root that is simple modulo $p_j$ lifts uniquely,
while a multiple root can lift to many roots or to none. The function
returns 0, without roots, if a list of roots would exceed $2^{32}$ entries
on 64-bit systems.

## Example {#example}

Modulo 7, $f = (x - 1)^3 (x - 2)(x^2 + 1)$ has the roots 1 and 2 (the
factor $x^2 + 1$ is irreducible since $7 \equiv 3 \bmod 4$). Modulo $120 =
2^3 \cdot 3 \cdot 5$, $x^2 - 1$ has $4 \cdot 2 \cdot 2 = 16$ roots, and
$x^2$ has the roots $a$ with $a^2 \equiv 0$, that is, the multiples of
$2^2 \cdot 3 \cdot 5 = 60$:

```c example name=nmod_poly_roots_composite
#include <stdlib.h>
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>

static int cmp(const void * a, const void * b)
{
    ulong x = *(const ulong *) a, y = *(const ulong *) b;
    return (x > y) - (x < y);
}

/* print the roots a (from the factors x - a) in increasing order,
   with the multiplicity if it is greater than 1 */
static void print_roots(const nmod_poly_factor_t r, ulong n)
{
    ulong a[64][2];
    slong i;
    for (i = 0; i < r->num; i++)
    {
        a[i][0] = (n - r->p[i].coeffs[0]) % n;
        a[i][1] = r->exp[i];
    }
    qsort(a, r->num, sizeof(a[0]), cmp);
    flint_printf("%wd roots mod %wu:", r->num, n);
    for (i = 0; i < r->num; i++)
    {
        flint_printf(" %wu", a[i][0]);
        if (a[i][1] > 1)
            flint_printf("^%wu", a[i][1]);
    }
    flint_printf("\n");
}

int main(void)
{
    nmod_poly_t f, g;
    nmod_poly_factor_t r;
    n_factor_t fac;
    slong i;

    nmod_poly_init(f, 7);
    nmod_poly_init(g, 7);
    nmod_poly_factor_init(r);

    nmod_poly_set_str(f, "2 7  6 1");          /* x - 1 */
    nmod_poly_pow(f, f, 3);
    nmod_poly_set_str(g, "2 7  5 1");          /* x - 2 */
    nmod_poly_mul(f, f, g);
    nmod_poly_set_str(g, "3 7  1 0 1");        /* x^2 + 1 */
    nmod_poly_mul(f, f, g);

    nmod_poly_roots(r, f, 1);
    for (i = 0; i < r->num; i++)
        flint_printf("root %wu with multiplicity %wd\n",
            nmod_neg(r->p[i].coeffs[0], f->mod), r->exp[i]);

    nmod_poly_clear(f);
    nmod_poly_init(f, 120);
    n_factor_init(&fac);
    n_factor(&fac, 120, 1);

    nmod_poly_set_str(f, "3 120  119 0 1");    /* x^2 - 1 */
    if (nmod_poly_roots_factored(r, f, 1, &fac))
        print_roots(r, 120);

    nmod_poly_set_str(f, "3 120  0 0 1");      /* x^2 */
    if (nmod_poly_roots_factored(r, f, 1, &fac))
        print_roots(r, 120);

    nmod_poly_clear(f);
    nmod_poly_clear(g);
    nmod_poly_factor_clear(r);
    return 0;
}
```
```output
root 2 with multiplicity 1
root 1 with multiplicity 3
16 roots mod 120: 1 11 19 29 31 41 49 59 61 71 79 89 91 101 109 119
2 roots mod 120: 0^2 60^2
```

The example reads only the constant coefficients of the factors, which
have the modulus of `f`. Both roots of $x^2$ have
multiplicity 2, since $(x - 60)^2 = x^2 - 120x +
3600 \equiv x^2 \pmod{120}$.
