---
title: Division
functions: fmpz_mpoly_divides fmpz_mpoly_divexact fmpz_mpoly_divrem fmpz_mpoly_div fmpz_mpoly_quasidivrem fmpz_mpoly_quasidiv fmpz_mpoly_divrem_ideal fmpz_mpoly_quasidivrem_ideal
---

Three kinds of division are provided for polynomials over $\Z$:

- *exact division*: `fmpz_mpoly_divides` decides whether $B$ divides $A$
  and returns the quotient if it does; `fmpz_mpoly_divexact` assumes that
  it does and throws otherwise;
- *division with remainder*: `fmpz_mpoly_divrem` and `fmpz_mpoly_div`
  compute $Q$ and $R$ with $A = QB + R$, where every coefficient division
  by $\operatorname{lc}(B)$ is a floor division;
- *quasi-division*: `fmpz_mpoly_quasidivrem` and `fmpz_mpoly_quasidiv`
  compute $Q$, $R$ and a positive integer $s$ with $sA = QB + R$ and no
  term of $R$ divisible by $\operatorname{lm}(B)$, which is division with
  remainder over $\Q$ with the denominators cleared.

The `_ideal` versions divide by several polynomials at once. All of them
use the heap algorithms of Monagan and Pearce, in which the terms of the
quotient are produced in order while a heap supplies the terms of
$A - QB$ that are still needed [@MonPea2007; @MonPea2011].

::: summary algorithms/multivariate/division :::

## Division with remainder over the integers {#divrem}

For univariate polynomials over a field, the division algorithm repeatedly
cancels the leading term of the dividend with a multiple of the divisor.
In several variables the same works with leading terms for a monomial
ordering: as long as the current dividend has a term $c\,m$ whose monomial
$m$ is divisible by $\operatorname{lm}(B)$, subtract $\frac{c}{\operatorname{lc}(B)}
\frac{m}{\operatorname{lm}(B)} B$; terms whose monomials are not
divisible go to the remainder. Since every step removes a term and only
adds smaller ones, this terminates, and the remainder has no monomial
divisible by $\operatorname{lm}(B)$ [@CLO2015, §2.3].

Over $\Z$ the coefficient division $c/\operatorname{lc}(B)$ need not be
exact. `fmpz_mpoly_divrem` performs the floor division $c = q
\operatorname{lc}(B) + r$, uses $q$ as the quotient coefficient and moves
the term $r\,m$ (if $r \ne 0$) to the remainder. The result satisfies $A =
QB + R$, and a term of $R$ whose monomial is divisible by
$\operatorname{lm}(B)$ has a coefficient of the same sign as
$\operatorname{lc}(B)$ and absolute value less than
$|\operatorname{lc}(B)|$. If $\operatorname{lc}(B) = \pm 1$, this is
ordinary division with remainder. If $B$ divides
$A$ exactly, $R = 0$ and $Q = A/B$ whatever $\operatorname{lc}(B)$ is,
and `fmpz_mpoly_div`, which does not form the remainder, is then the
fastest way to compute the quotient.

`fmpz_mpoly_quasidivrem` instead multiplies the dividend by just enough
to make each coefficient division exact: whenever a quotient coefficient
$c/\operatorname{lc}(B)$ is not an integer, everything computed so far is
scaled by $|\operatorname{lc}(B)|/\gcd(c, \operatorname{lc}(B))$. The
result satisfies $sA = QB + R$ with no monomial of $R$ divisible by
$\operatorname{lm}(B)$, and $s$ divides a power of
$\operatorname{lc}(B)$. Dividing by $s$ gives the quotient and remainder
of the division in $\Q[x_0, \ldots, x_{n-1}]$.

With several divisors $B_0, \ldots, B_{k-1}$ (`fmpz_mpoly_divrem_ideal`,
`fmpz_mpoly_quasidivrem_ideal`), each step uses the first $B_i$ whose
leading monomial divides the current term. The remainder then depends on
the order of the divisors unless they form a Gröbner basis.

## Exact division {#exact}

`fmpz_mpoly_divides` runs the division algorithm but stops as soon as a
term appears that would go to the remainder: a leading term not divisible
by $\operatorname{lm}(B)$ or a coefficient not divisible by
$\operatorname{lc}(B)$. Before starting, it returns 0 if $B$ has larger
degree than $A$ in some variable. Since an exact quotient has no larger
exponents than $A$, the computation is done with the packing width of
$A$, and an exponent overflow also proves that the division is not exact.
In the common case where $B$ does not divide $A$ it therefore often
returns after a few terms. It throws `FLINT_DIVZERO` for a zero divisor.

If FLINT knows that the machine has a strongly ordered memory model
(`FLINT_KNOW_STRONG_ORDER`, true on x86), dividends of more than 1024
terms are divided with threads, following [@MonPea2010]: the monomials of
$A$ are cut into ranges at boundaries chosen by
`mpoly_divides_select_exps`; each worker subtracts the products of the
quotient terms found so far from its range and, when its range becomes the
leading one, produces the next quotient terms.

The array methods `fmpz_mpoly_divides_array` and `fmpz_mpoly_divrem_array`
are not used by these functions; they can be faster for dense inputs with
small exponents.

## Examples {#examples}

Dividing $A = 3x^2 + 5xy + 7$ by $B = 2x + y$ in lex order with $x > y$:
the floor division of 3 by 2 leaves $x^2$ in the remainder, while the
quasi-division scales $A$ by $s = 4$ and leaves no multiple of $x$ in the
remainder. Divisibility is then tested on $AB$ and $AB + 1$.

```c example name=fmpz_mpoly_division
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>

int main(void)
{
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t A, B, Q, R, T;
    fmpz_t s;
    const char * vars[] = {"x", "y"};

    fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);
    fmpz_mpoly_init(A, ctx); fmpz_mpoly_init(B, ctx);
    fmpz_mpoly_init(Q, ctx); fmpz_mpoly_init(R, ctx); fmpz_mpoly_init(T, ctx);
    fmpz_init(s);

    fmpz_mpoly_set_str_pretty(A, "3*x^2 + 5*x*y + 7", vars, ctx);
    fmpz_mpoly_set_str_pretty(B, "2*x + y", vars, ctx);

    /* A = Q*B + R with floor division of the coefficients by lc(B) = 2 */
    fmpz_mpoly_divrem(Q, R, A, B, ctx);
    flint_printf("divrem:      Q = "); fmpz_mpoly_print_pretty(Q, vars, ctx);
    flint_printf(", R = "); fmpz_mpoly_print_pretty(R, vars, ctx); flint_printf("\n");

    /* s*A = Q*B + R with no term of R divisible by x */
    fmpz_mpoly_quasidivrem(s, Q, R, A, B, ctx);
    flint_printf("quasidivrem: s = "); fmpz_print(s);
    flint_printf(", Q = "); fmpz_mpoly_print_pretty(Q, vars, ctx);
    flint_printf(", R = "); fmpz_mpoly_print_pretty(R, vars, ctx); flint_printf("\n");
    fmpz_mpoly_mul(T, Q, B, ctx);
    fmpz_mpoly_add(T, T, R, ctx);
    fmpz_mpoly_scalar_mul_fmpz(Q, A, s, ctx);
    flint_printf("check: %d\n", fmpz_mpoly_equal(T, Q, ctx));

    /* exact division */
    fmpz_mpoly_mul(T, A, B, ctx);
    flint_printf("divides(A*B, B): %d", fmpz_mpoly_divides(Q, T, B, ctx));
    flint_printf(", quotient == A: %d\n", fmpz_mpoly_equal(Q, A, ctx));
    fmpz_mpoly_add_ui(T, T, 1, ctx);
    flint_printf("divides(A*B + 1, B): %d\n", fmpz_mpoly_divides(Q, T, B, ctx));

    fmpz_clear(s);
    fmpz_mpoly_clear(A, ctx); fmpz_mpoly_clear(B, ctx);
    fmpz_mpoly_clear(Q, ctx); fmpz_mpoly_clear(R, ctx); fmpz_mpoly_clear(T, ctx);
    fmpz_mpoly_ctx_clear(ctx);
    return 0;
}
```
```output
divrem:      Q = x+2*y, R = x^2-2*y^2+7
quasidivrem: s = 4, Q = 6*x+7*y, R = -7*y^2+28
check: 1
divides(A*B, B): 1, quotient == A: 1
divides(A*B + 1, B): 0
```

In `divrem`, the first quotient term is $\lfloor 3/2 \rfloor x = x$, which
leaves $x^2 + 4xy + 7$ after subtracting $xB$; the term $x^2$ has
coefficient $1 < 2$ and goes to the remainder, and $4xy$ gives the
quotient term $2y$. In `quasidivrem`, both the term $3x^2$ and the later
term $7xy$ force a factor 2, so $s = 4 = \operatorname{lc}(B)^2$.
