---
title: Representation and arithmetic of rational functions
summary: The canonical form of an fmpz_poly_q_t, how Henrici's algorithms produce canonical sums and products with small gcds, and how to build a rational function from its numerator and denominator.
order: 1
---

An `fmpz_poly_q_t` represents $f \in \mathbb{Q}(x)$ as a quotient $a/b$ of
integer polynomials with $\gcd(a, b) = 1$ in $\Z[x]$ and $b$ having
positive leading coefficient. Any $f$ has exactly one such representation:
clearing denominators of a representation over $\mathbb{Q}[x]$ gives one over
$\Z[x]$, and two coprime representations $a/b = a'/b'$ satisfy $a b' = a'
b$, so $b \mid b'$ and $b' \mid b$ by Gauss's lemma, hence $b' = \pm b$.
Since the gcd in $\Z[x]$ includes the content, $(2x + 2)/(4x)$ is not
canonical, while $(x + 1)/(2x)$ is.

Canonical forms make equality testing a comparison of coefficients, but
reducing a quotient costs a polynomial gcd in $\Z[x]$, which is much more
expensive than a product. The arithmetic functions therefore never
reduce a full product. They use Henrici's algorithms [@Hen1956], which
exploit the coprimality of the inputs to cancel common factors before
multiplying.

## Henrici's algorithms {#henrici}

Let $a/b$ and $c/d$ be canonical.

::: lemma Products {#lem-product}
With $g_1 = \gcd(a, d)$ and $g_2 = \gcd(c, b)$,
$$
\frac{a}{b} \cdot \frac{c}{d} = \frac{(a/g_1)(c/g_2)}{(b/g_2)(d/g_1)},
$$
and the right-hand side is canonical up to the sign of the denominator.
:::

::: proof
An irreducible $\pi$ dividing both the numerator and the denominator
divides one of $a/g_1$, $c/g_2$ and one of $b/g_2$, $d/g_1$. It cannot
divide both $a/g_1$ and $b/g_2$, since $\gcd(a, b) = 1$, and it cannot
divide both $a/g_1$ and $d/g_1$, by the definition of $g_1$; the other two
cases are symmetric.
:::

::: lemma Sums {#lem-sum}
With $g = \gcd(b, d)$, $b = g b'$, $d = g d'$ and $t = a d' + c b'$,
$$
\frac{a}{b} + \frac{c}{d} = \frac{t}{b' d' g}, \qquad \gcd(t, b' d' g) =
\gcd(t, g).
$$
:::

::: proof
$\gcd(t, b') = \gcd(a d', b') = 1$, because $\gcd(a, b') = 1$ and
$\gcd(d', b') = 1$; likewise $\gcd(t, d') = 1$. So every common factor of
$t$ and $b' d' g$ divides $g$.
:::

So a sum costs one gcd of the denominators and one gcd of $t$ with the
(usually small) $g$, and a product two gcds of polynomials of the size of
the inputs, instead of a gcd of the full numerator and denominator. When
the denominators are coprime, which is the common case, no cancellation
happens at all. `fmpz_poly_q_add`, `fmpz_poly_q_sub`, `fmpz_poly_q_mul`
and `fmpz_poly_q_div` implement these formulas, with shortcuts for zero
operands, for denominators equal to 1 and for constant denominators,
where the arithmetic of `fmpq_poly` is used. The derivative is computed
with a similar saving: $(a/b)' = (a' b_1 - a b'/d)/(b b_1)$ with $d =
\gcd(b, b')$ and $b_1 = b/d$, whose numerator and denominator can only
share a constant factor [@Hor1972].

## Working with numerators and denominators {#numref}

`fmpz_poly_q_numref` and `fmpz_poly_q_denref` give access to the two
`fmpz_poly_struct`s, so any `fmpz_poly` function can be used to build a
rational function, after which `fmpz_poly_q_canonicalise` restores the
canonical form. Strings in the format of `fmpz_poly_q_set_str` are an
alternative for small inputs.

```c example name=fmpz_poly_q_basic
#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_q.h>

int main(void)
{
    fmpz_poly_q_t f, g, h;
    fmpq_t a, v;

    fmpz_poly_q_init(f);
    fmpz_poly_q_init(g);
    fmpz_poly_q_init(h);
    fmpq_init(a);
    fmpq_init(v);

    /* f = (2x^2 - 2)/(4x^2 + 4x), built directly, then canonicalised */
    fmpz_poly_set_str(fmpz_poly_q_numref(f), "3  -2 0 2");
    fmpz_poly_set_str(fmpz_poly_q_denref(f), "3  0 4 4");
    fmpz_poly_q_canonicalise(f);
    flint_printf("f = "); fmpz_poly_q_print_pretty(f, "x");
    flint_printf("  ("); fmpz_poly_q_print(f); flint_printf(")\n");

    fmpz_poly_q_set_str(g, "2  1 3/1  2");          /* (3x + 1)/2 */
    fmpz_poly_q_mul(h, f, g);
    flint_printf("f g = "); fmpz_poly_q_print_pretty(h, "x"); flint_printf("\n");
    fmpz_poly_q_add(h, f, g);
    flint_printf("f + g = "); fmpz_poly_q_print_pretty(h, "x"); flint_printf("\n");
    fmpz_poly_q_derivative(h, f);
    flint_printf("f' = "); fmpz_poly_q_print_pretty(h, "x"); flint_printf("\n");

    fmpq_set_si(a, 1, 3);
    if (fmpz_poly_q_evaluate_fmpq(v, f, a) == 0)
    {
        flint_printf("f(1/3) = "); fmpq_print(v); flint_printf("\n");
    }
    fmpq_set_si(a, 0, 1);
    flint_printf("f(0) defined: %d\n", fmpz_poly_q_evaluate_fmpq(v, f, a) == 0);

    fmpz_poly_q_clear(f);
    fmpz_poly_q_clear(g);
    fmpz_poly_q_clear(h);
    fmpq_clear(a);
    fmpq_clear(v);
    return 0;
}
```
```output
f = (x-1)/(2*x)  (2  -1 1/2  0 2)
f g = (3*x^2-2*x-1)/(4*x)
f + g = (3*x^2+2*x-1)/(2*x)
f' = 1/(2*x^2)
f(1/3) = -1
f(0) defined: 0
```

The common factor $2(x + 1)$ of the numerator and denominator of $f$ is
removed by `fmpz_poly_q_canonicalise`. In $f + g$ the denominators $2x$ and
$2$ have the gcd $g = 2$, and since the numerator $3x^2 + 2x - 1$ is
coprime to 2 in $\Z[x]$, nothing more cancels.
