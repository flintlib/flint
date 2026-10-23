---
title: Resultants and discriminants
functions: nmod_poly_resultant nmod_poly_resultant_euclidean nmod_poly_resultant_hgcd nmod_poly_discriminant
---

The resultant $\operatorname{res}(f, g)$ is the determinant of the
Sylvester matrix of $f$ and $g$; it vanishes exactly when $f$ and $g$
have a common factor, and for $f = a \prod (x - \alpha_i)$ it equals $a^{\deg
g} \prod_i g(\alpha_i)$. The discriminant of $f$ is, up to sign and a power
of the leading coefficient, the resultant of $f$ and $f'$, and it vanishes
exactly when $f$ has a repeated root. Resultants over $\Z/p\Z$ for many
primes $p$ are the usual way to compute resultants of integer polynomials
by Chinese remaindering. These functions require a prime modulus.

## Algorithms {#algorithms}

The resultant is computed along the Euclidean remainder sequence: for $r =
f \bmod g$,
$$
\operatorname{res}(f, g) = (-1)^{\deg f \deg g}\,
\operatorname{lc}(g)^{\deg f - \deg r} \operatorname{res}(g, r),
$$
and $\operatorname{res}(f, c) = c^{\deg f}$ for a constant $c$
[@GG2003, §6]. `nmod_poly_resultant_euclidean` follows this recursion in
$O(\deg f \deg g)$ operations. `nmod_poly_resultant_hgcd` performs the same
divisions with the [half-GCD algorithm](algorithms/polynomials/half-gcd)
in $O(M(n) \log n)$ operations; because a recursive half-GCD call only
knows the top coefficients of the remainders, the factors contributed by
its last division are recorded and applied after the remainders have been
completed. `nmod_poly_resultant` handles small degrees by explicit
formulas, and otherwise chooses between the two at the same cutoff as
`nmod_poly_gcd`, using Montgomery arithmetic for the Euclidean algorithm
under the same conditions.

The discriminant is computed as
$$
\operatorname{disc}(f) = (-1)^{d(d-1)/2}\, \operatorname{lc}(f)^{d -
\deg f' - 2} \operatorname{res}(f, f'),
$$
which is $(-1)^{d(d-1)/2} \operatorname{res}(f, f') / \operatorname{lc}(f)$
when $\deg f' = d - 1$, that is, when the characteristic does not divide
$d$.

## Examples {#examples}

For the cubic $x^3 + ax + b$ the discriminant is $-4a^3 - 27b^2$. Modulo
$p = 101$, the program checks this formula, finds the values of $b$ for
which $x^3 + 3x + b$ has a repeated root, and compares the two resultant
algorithms on long random input:

```c example name=nmod_poly_resultant_disc
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>

int main(void)
{
    nmod_poly_t f, g, h;
    flint_rand_t state;
    nmod_t mod;
    ulong a = 3, b, d;

    nmod_init(&mod, 101);
    flint_rand_init(state);
    nmod_poly_init_mod(f, mod);
    nmod_poly_init_mod(g, mod);
    nmod_poly_init_mod(h, mod);

    nmod_poly_set_coeff_ui(f, 3, 1);
    nmod_poly_set_coeff_ui(f, 1, a);
    nmod_poly_set_coeff_ui(f, 0, 5);
    d = nmod_sub(nmod_neg(nmod_mul(4, nmod_pow_ui(a, 3, mod), mod), mod),
                 nmod_mul(27, 25, mod), mod);
    flint_printf("disc = %wu, -4a^3 - 27b^2 = %wu\n",
                 nmod_poly_discriminant(f), d);

    flint_printf("repeated root for b =");
    for (b = 0; b < mod.n; b++)
    {
        nmod_poly_set_coeff_ui(f, 0, b);
        if (nmod_poly_discriminant(f) == 0)
            flint_printf(" %wu", b);
    }
    flint_printf("\n");

    /* res(f, g) = 0 exactly when gcd(f, g) != 1 */
    nmod_poly_set_coeff_ui(f, 0, 5);
    nmod_poly_set_str(g, "2 101  1 1");               /* x + 1 */
    flint_printf("res(f, x + 1) = %wu = -f(-1) = %wu\n",
                 nmod_poly_resultant(f, g),
                 nmod_neg(nmod_poly_evaluate_nmod(f, mod.n - 1), mod));

    nmod_poly_rand(g, state, 3000);
    nmod_poly_rand(h, state, 2500);
    flint_printf("euclidean == hgcd: %d\n",
                 nmod_poly_resultant_euclidean(g, h) ==
                 nmod_poly_resultant_hgcd(g, h));

    nmod_poly_clear(f);
    nmod_poly_clear(g);
    nmod_poly_clear(h);
    flint_rand_clear(state);
    return 0;
}
```
```output
disc = 25, -4a^3 - 27b^2 = 25
repeated root for b = 20 81
res(f, x + 1) = 100 = -f(-1) = 100
euclidean == hgcd: 1
```
