---
title: Getting started with multivariate polynomials
summary: A short program that creates the ring Q[x, y], parses polynomials, tests divisibility and divides with remainder, and notes on monomial orderings and coefficient rings that are not fields.
order: 1
---

A `gr_mpoly_ctx_t` is a generic ring: it is created from a coefficient
ring, a number of variables and a monomial ordering, and its elements
can be manipulated with the functions of `gr.h` (`gr_set_str`, `gr_mul`,
`gr_println`, ...) as well as with the `gr_mpoly` functions, which take
the same context. The program below works in $\mathbb{Q}[x, y]$ with the
degree-lexicographic ordering.

```c example name=gr_mpoly_division
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_mpoly.h>

int main(void)
{
    gr_ctx_t QQ;
    gr_mpoly_ctx_t R;
    gr_mpoly_t A, B, Q, Rem;
    const char * vars[] = { "x", "y" };
    int status;

    gr_ctx_init_fmpq(QQ);
    gr_mpoly_ctx_init(R, QQ, 2, ORD_DEGLEX);       /* Q[x, y] */
    GR_MUST_SUCCEED(gr_ctx_set_gen_names(R, vars));

    gr_mpoly_init(A, R);
    gr_mpoly_init(B, R);
    gr_mpoly_init(Q, R);
    gr_mpoly_init(Rem, R);

    GR_MUST_SUCCEED(gr_set_str(A, "(x + y)^3 - 1", R));
    GR_MUST_SUCCEED(gr_set_str(B, "x + y - 1", R));
    flint_printf("A = "); gr_println(A, R);

    status = gr_mpoly_divides(Q, A, B, R);
    flint_printf("divides: status %d, Q = ", status); gr_println(Q, R);

    GR_MUST_SUCCEED(gr_set_str(B, "x*y + 1", R));
    GR_MUST_SUCCEED(gr_mpoly_divrem(Q, Rem, A, B, R));
    flint_printf("A = Q*(x*y + 1) + R with Q = "); gr_print(Q, R);
    flint_printf(", R = "); gr_println(Rem, R);

    status = gr_mpoly_divides(Q, A, B, R);
    flint_printf("divides: status %d\n", status);

    gr_mpoly_clear(A, R);
    gr_mpoly_clear(B, R);
    gr_mpoly_clear(Q, R);
    gr_mpoly_clear(Rem, R);
    gr_ctx_clear(R);
    gr_ctx_clear(QQ);
    return 0;
}
```
```output
A = x^3 + 3*x^2*y + 3*x*y^2 + y^3 - 1
divides: status 0, Q = x^2 + 2*x*y + y^2 + x + y + 1
A = Q*(x*y + 1) + R with Q = 3*x + 3*y, R = x^3 + y^3 - 3*x - 3*y - 1
divides: status 1
```

The generator names default to `x1`, `x2`, ...; `gr_ctx_set_gen_names`
changes them, and the parser then accepts expressions in the new names.
`gr_mpoly_divides` decides divisibility and returns `GR_DOMAIN` (1)
when $B$ does not divide $A$, so that a status code, not a sentinel
value of $Q$, carries the answer.

Division with remainder depends on the monomial ordering: the leading
monomial of `x*y + 1` is $xy$ for every ordering here, and the
remainder is reduced in the sense that none of its terms is divisible
by $xy$. With several divisors (`gr_mpoly_divrem_ideal`) the remainder
also depends on the order of the divisors unless they form a Gröbner
basis [@CLO2015].

Over a coefficient ring that is not a field, such as $\Z$, the
leading coefficient of the divisor need not be invertible.
`gr_mpoly_divrem` then returns `GR_DOMAIN` unless every division of
coefficients that occurs is exact, while the *weak* division functions
(`gr_mpoly_divrem_weak`) divide such coefficients with remainder in
the coefficient ring and keep the remainder, so that $A = QB + R$ still
holds but terms of $R$ may be divisible by the leading monomial of
$B$. The functions and their
algorithms, including the heap algorithms of Monagan and Pearce
[@MonPea2007; @MonPea2009] used for multiplication and division, are
described in [](gr_mpoly).
