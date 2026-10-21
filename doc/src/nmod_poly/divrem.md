---
title: Division with remainder
functions: nmod_poly_divrem nmod_poly_div nmod_poly_rem nmod_poly_divexact nmod_poly_divides
---

For polynomials $A$ and $B \ne 0$ whose leading coefficient is invertible
modulo $n$, there are unique $Q$ and $R$ with $A = BQ + R$ and $\deg R <
\deg B$. `nmod_poly_divrem` computes both; `nmod_poly_div` and
`nmod_poly_rem` compute one of them, `nmod_poly_divexact` computes $Q$
when $R$ is known to be zero, and `nmod_poly_divides` decides whether $R =
0$. The modulus need not be prime, but if the leading coefficient of $B$
is not a unit, an exception is thrown even when the division would be
possible.

## Algorithms {#algorithms}

Let $\mathit{lenQ} = \deg A - \deg B + 1$ be the length of the quotient.

- *Schoolbook division* computes the quotient coefficients from the top,
  subtracting a multiple of $B$ for each, with $O(\mathit{lenQ} \cdot
  \mathit{lenB})$ operations. FLINT keeps the partial remainder unreduced
  in one to three words per coefficient and reduces a coefficient only
  when it becomes the leading one, so the inner loop is a plain
  multiply-accumulate.
- *Newton division* uses the identity $\operatorname{rev}(Q) =
  \operatorname{rev}(A) / \operatorname{rev}(B) \bmod x^{\mathit{lenQ}}$
  for the reversed polynomials (see [](nmod_poly/modular-arithmetic)):
  the quotient is a power series quotient, computed by
  [Newton iteration](algorithms/polynomials/newton) in a constant number
  of multiplications of length $\mathit{lenQ}$, and the remainder is $A -
  BQ$, of which only the low $\mathit{lenB} - 1$ coefficients are needed.

`nmod_poly_divrem` uses schoolbook division when the quotient has at most
80 coefficients or the divisor is short (fewer than 400 coefficients for
moduli of at most 27 bits, 80 otherwise, in builds with AVX2), and Newton
division otherwise. `nmod_poly_div` switches to Newton division earlier,
when both lengths exceed 15, since it does not need the remainder.
`nmod_poly_rem` handles the case $\mathit{lenQ} = 2$, the typical step of
the Euclidean algorithm, with a dedicated routine. Division by a sparse
modulus is fastest with `nmod_poly_divrem_newton_n_preinv`, which reduces
directly when the divisor has at most five nonzero lower terms.

`nmod_poly_divexact` exploits exactness: the quotient is determined both
by the top coefficients of $A$ (division from the top) and by its bottom
coefficients (power series division by $B$, if $B(0)$ is a unit), and for
medium and long inputs it computes half of the quotient from each end,
which roughly halves the work. `nmod_poly_divides` computes the quotient
with `nmod_poly_div` and then checks $A - BQ$ only in its low coefficients,
stopping at the first nonzero one.

For repeated division by the same $B$, precompute the inverse of the
reversed divisor and use `nmod_poly_divrem_newton_n_preinv`; see
[](nmod_poly/modular-arithmetic).

## Examples {#examples}

Modulo 13, the program builds $A = (x^3 + 2x + 5)(x^2 + 3) + x + 1$ and
recovers the quotient and remainder, divides exactly, tests divisibility
and removes a repeated factor:

```c example name=nmod_poly_divrem_basic
#include <flint/flint.h>
#include <flint/nmod_poly.h>

int main(void)
{
    nmod_poly_t a, b, c, q, r;

    nmod_poly_init(a, 13);
    nmod_poly_init(b, 13);
    nmod_poly_init(c, 13);
    nmod_poly_init(q, 13);
    nmod_poly_init(r, 13);

    nmod_poly_set_str(c, "4 13  5 2 0 1");    /* x^3 + 2x + 5 */
    nmod_poly_set_str(b, "3 13  3 0 1");      /* x^2 + 3 */
    nmod_poly_mul(a, c, b);
    nmod_poly_add_ui(a, a, 1);
    nmod_poly_set_coeff_ui(a, 1, nmod_poly_get_coeff_ui(a, 1) + 1);

    nmod_poly_divrem(q, r, a, b);
    flint_printf("q = "); nmod_poly_print_pretty(q, "x");
    flint_printf(", r = "); nmod_poly_print_pretty(r, "x"); flint_printf("\n");

    nmod_poly_sub(a, a, r);                   /* now b divides a */
    nmod_poly_divexact(q, a, b);
    flint_printf("divexact: "); nmod_poly_print_pretty(q, "x"); flint_printf("\n");

    flint_printf("b | a: %d, ", nmod_poly_divides(q, a, b));
    nmod_poly_add_ui(a, a, 1);
    flint_printf("b | a + 1: %d\n", nmod_poly_divides(q, a, b));

    /* a = b^3 c */
    nmod_poly_pow(a, b, 3);
    nmod_poly_mul(a, a, c);
    flint_printf("multiplicity of b: %wu, cofactor ", nmod_poly_remove(a, b));
    nmod_poly_print_pretty(a, "x"); flint_printf("\n");

    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(c);
    nmod_poly_clear(q);
    nmod_poly_clear(r);
    return 0;
}
```
```output
q = x^3+2*x+5, r = x+1
divexact: x^3+2*x+5
b | a: 1, b | a + 1: 0
multiplicity of b: 3, cofactor x^3+2*x+5
```

Division by a non-monic divisor modulo a composite number works as long
as its leading coefficient is a unit. Modulo 15, the program divides by
$7x^2 + 1$ ($7 \cdot 13 \equiv 1$) and checks $A = BQ + R$:

```c example name=nmod_poly_divrem_composite
#include <flint/flint.h>
#include <flint/nmod_poly.h>

int main(void)
{
    nmod_poly_t a, b, q, r, t;

    nmod_poly_init(a, 15);
    nmod_poly_init(b, 15);
    nmod_poly_init(q, 15);
    nmod_poly_init(r, 15);
    nmod_poly_init(t, 15);

    nmod_poly_set_str(a, "6 15  1 2 3 4 5 6");
    nmod_poly_set_str(b, "3 15  1 0 7");

    nmod_poly_divrem(q, r, a, b);
    flint_printf("q = "); nmod_poly_print_pretty(q, "x");
    flint_printf(", r = "); nmod_poly_print_pretty(r, "x"); flint_printf("\n");

    nmod_poly_mul(t, b, q);
    nmod_poly_add(t, t, r);
    flint_printf("b q + r == a: %d\n", nmod_poly_equal(t, a));

    nmod_poly_clear(a);
    nmod_poly_clear(b);
    nmod_poly_clear(q);
    nmod_poly_clear(r);
    nmod_poly_clear(t);
    return 0;
}
```
```output
q = 3*x^3+5*x^2+13*x+4, r = 4*x+12
b q + r == a: 1
```
