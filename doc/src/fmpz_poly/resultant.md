---
title: Resultants and discriminants
functions: fmpz_poly_resultant fmpz_poly_resultant_modular fmpz_poly_resultant_modular_div fmpz_poly_resultant_euclidean fmpz_poly_discriminant
---

The resultant of $f = a \prod_{i=1}^m (x - \alpha_i)$ and $g = b
\prod_{j=1}^n (x - \beta_j)$ is

$$
\operatorname{res}(f, g) = a^n b^m \prod_{i,j} (\alpha_i - \beta_j)
= a^n \prod_{i=1}^m g(\alpha_i) = (-1)^{mn} b^m \prod_{j=1}^n f(\beta_j),
$$

the determinant of the $(m + n) \times (m + n)$ Sylvester matrix of $f$
and $g$. It is an integer that vanishes exactly when $f$ and $g$ have a
common root, and it satisfies $\operatorname{res}(f, g) = (-1)^{mn}
\operatorname{res}(g, f)$, $\operatorname{res}(f, gh) =
\operatorname{res}(f, g) \operatorname{res}(f, h)$ and
$\operatorname{res}(cf, g) = c^n \operatorname{res}(f, g)$. The
discriminant $\operatorname{disc}(f) = (-1)^{m(m-1)/2}
\operatorname{res}(f, f') / a$ vanishes exactly when $f$ has a
repeated root. FLINT defines the resultant to be zero if one of the
polynomials is zero.

::: summary algorithms/polynomials/resultant-integers :::

## Algorithms {#algorithms}

`fmpz_poly_resultant` chooses between two algorithms. For short
polynomials with small coefficients, precisely when the shorter has
length $\ell_2 \le 144$ and $\ell_2^3 (b_1 + b_2) \le 6 \cdot 10^6$ for
coefficient sizes $b_1, b_2$ bits, it runs the subresultant algorithm
(`fmpz_poly_resultant_euclidean`), a pseudo-remainder sequence whose
coefficients are kept small by exact divisions. Otherwise it uses the
multimodular algorithm (`fmpz_poly_resultant_modular`): Hadamard's
inequality bounds the resultant by $\|f\|_2^{n} \|g\|_2^{m}$, and the
resultant is computed modulo enough primes of 64 bits with
`nmod_poly_resultant` (half-GCD based for long polynomials) and
reconstructed by the Chinese remainder theorem. The number of primes
is fixed in advance by the bound, which is usually close to the actual
size of the resultant.

When the size of the result is known for other reasons, or a large
known factor can be removed, `fmpz_poly_resultant_modular_div`
computes $\operatorname{res}(f, g)/d$ for a given exact divisor $d$
with only as many primes as a given bit bound for the quotient
requires; the bound should be one bit more than the size of the
quotient. The discriminant is computed as a resultant with the
derivative.

## Examples {#examples}

The resultant of $x^2 - 2$ and $x^2 - 3$ is $\prod (\pm\sqrt 2 \mp
\sqrt 3)$, a product over the four pairs of roots, and the discriminant
of $x^3 - x - 1$ is $-23$, the discriminant of the cubic field it
defines. The resultant changes sign when two polynomials of odd degree
are swapped.

```c example name=fmpz_poly_resultant_small
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

int main(void)
{
    fmpz_poly_t f, g;
    fmpz_t r;

    fmpz_poly_init(f);
    fmpz_poly_init(g);
    fmpz_init(r);

    fmpz_poly_set_str(f, "3  -2 0 1");        /* x^2 - 2 */
    fmpz_poly_set_str(g, "3  -3 0 1");        /* x^2 - 3 */
    fmpz_poly_resultant(r, f, g);
    flint_printf("res(x^2 - 2, x^2 - 3) = "); fmpz_print(r); flint_printf("\n");

    fmpz_poly_set_str(f, "4  -1 -1 0 1");     /* x^3 - x - 1 */
    fmpz_poly_discriminant(r, f);
    flint_printf("disc(x^3 - x - 1) = "); fmpz_print(r); flint_printf("\n");

    fmpz_poly_set_str(g, "2  5 2");           /* 2x + 5 */
    fmpz_poly_resultant(r, f, g);
    flint_printf("res(f, 2x + 5) = "); fmpz_print(r);
    fmpz_poly_resultant(r, g, f);
    flint_printf(", res(2x + 5, f) = "); fmpz_print(r); flint_printf("\n");

    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    fmpz_clear(r);
    return 0;
}
```
```output
res(x^2 - 2, x^2 - 3) = 1
disc(x^3 - x - 1) = -23
res(f, 2x + 5) = 113, res(2x + 5, f) = -113
```

Indeed $(\sqrt 2 - \sqrt 3)(\sqrt 2 + \sqrt 3)$ times its conjugate
under $\sqrt 2 \mapsto -\sqrt 2$ is $(-1)(-1) = 1$, and
$\operatorname{res}(2x + 5, f) = 2^3 f(-5/2) = -113$.

For larger inputs the two algorithms give the same value. The second
part divides out a known factor: $\operatorname{res}(f, gh) =
\operatorname{res}(f, g)\operatorname{res}(f, h)$, so knowing
$\operatorname{res}(f, h)$ and a bound for the quotient saves primes.

```c example name=fmpz_poly_resultant_large
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

int main(void)
{
    flint_rand_t state;
    fmpz_poly_t f, g, h, gh;
    fmpz_t r1, r2, d, q;

    flint_rand_init(state);
    fmpz_poly_init(f);
    fmpz_poly_init(g);
    fmpz_poly_init(h);
    fmpz_poly_init(gh);
    fmpz_init(r1);
    fmpz_init(r2);
    fmpz_init(d);
    fmpz_init(q);

    fmpz_poly_randtest(f, state, 60, 30);
    fmpz_poly_randtest(g, state, 50, 30);
    fmpz_poly_set_coeff_ui(f, 59, 1);
    fmpz_poly_set_coeff_ui(g, 49, 1);

    fmpz_poly_resultant_modular(r1, f, g);
    fmpz_poly_resultant_euclidean(r2, f, g);
    flint_printf("modular == euclidean: %d, %wd bits\n",
                 fmpz_equal(r1, r2), fmpz_bits(r1));

    fmpz_poly_set_str(h, "3  1 1 1");         /* x^2 + x + 1 */
    fmpz_poly_mul(gh, g, h);
    fmpz_poly_resultant(d, f, h);
    fmpz_poly_resultant_modular_div(q, f, gh, d, fmpz_bits(r1) + 1);
    flint_printf("res(f, gh) / res(f, h) == res(f, g): %d\n", fmpz_equal(q, r1));

    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    fmpz_poly_clear(h);
    fmpz_poly_clear(gh);
    fmpz_clear(r1);
    fmpz_clear(r2);
    fmpz_clear(d);
    fmpz_clear(q);
    flint_rand_clear(state);
    return 0;
}
```
```output
modular == euclidean: 1, 2940 bits
res(f, gh) / res(f, h) == res(f, g): 1
```
