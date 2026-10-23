---
title: Large algebraic expressions
summary: Parsing expressions with radicals in the generic rings interface and deciding equalities with qqbar and ca, and how the two representations divide the work differently on an expression with thousands of operations; the example program huge_expr.
order: 4
---

Symbolic computation sometimes produces algebraic numbers as large
expressions with nested radicals, and deciding whether two of them are
equal is a common and difficult task for computer algebra systems. FLINT
can evaluate such expressions either as canonical algebraic numbers
([](qqbar)) or as elements of formal fields ([](ca)).

## Parsing and comparing {#parsing}

The generic rings interface parses expressions with `gr_set_str`, which
understands integers, the arithmetic operators, `^` with rational
exponents and functions such as `sqrt`, and evaluates them in the given
context. The program below checks two denesting identities, Ramanujan's

$$
\sqrt[3]{\sqrt[3]{2} - 1} = \sqrt[3]{1/9} - \sqrt[3]{2/9} + \sqrt[3]{4/9}
$$

and $\sqrt{2 + \sqrt 3} = (\sqrt 6 + \sqrt 2)/2$, in the field of
algebraic numbers with `qqbar_t` elements (`gr_ctx_init_complex_qqbar`)
and with `ca_t` elements (`gr_ctx_init_complex_algebraic_ca`). With
[](qqbar), each side is computed with its minimal polynomial and the
comparison is immediate; with [](ca), each side is an element of a field
with several generators, and `gr_equal` calls the zero test, which
evaluates the difference as an algebraic number (see
[](algorithms/exact/zero-testing)).

```c example name=tutorial_exact_expressions
#include <flint/calcium.h>
#include <flint/gr.h>

static void
check(const char * lhs, const char * rhs, gr_ctx_t ctx)
{
    gr_ptr a, b;
    int status;

    GR_TMP_INIT2(a, b, ctx);
    status = gr_set_str(a, lhs, ctx);
    status |= gr_set_str(b, rhs, ctx);
    if (status == GR_SUCCESS)
    {
        flint_printf("%s = %s: ", lhs, rhs);
        truth_print(gr_equal(a, b, ctx));
        flint_printf("\n");
    }
    else
        flint_printf("evaluation failed\n");
    GR_TMP_CLEAR2(a, b, ctx);
}

int main(void)
{
    gr_ctx_t QQbar, CA;

    gr_ctx_init_complex_qqbar(QQbar);
    gr_ctx_init_complex_algebraic_ca(CA);

    check("(2^(1/3)-1)^(1/3)", "(1/9)^(1/3)-(2/9)^(1/3)+(4/9)^(1/3)", QQbar);
    check("(2^(1/3)-1)^(1/3)", "(1/9)^(1/3)-(2/9)^(1/3)+(4/9)^(1/3)", CA);
    check("sqrt(2+sqrt(3))", "(sqrt(6)+sqrt(2))/2", QQbar);
    check("sqrt(2+sqrt(3))", "(sqrt(6)+sqrt(2))/2", CA);

    gr_ctx_clear(QQbar);
    gr_ctx_clear(CA);
    return 0;
}
```

```output
(2^(1/3)-1)^(1/3) = (1/9)^(1/3)-(2/9)^(1/3)+(4/9)^(1/3): T_TRUE
(2^(1/3)-1)^(1/3) = (1/9)^(1/3)-(2/9)^(1/3)+(4/9)^(1/3): T_TRUE
sqrt(2+sqrt(3)) = (sqrt(6)+sqrt(2))/2: T_TRUE
sqrt(2+sqrt(3)) = (sqrt(6)+sqrt(2))/2: T_TRUE
```

## A huge expression {#huge}

The example program `huge_expr` solves a problem posted on the Sage
help forum: two algebraic numbers $N$ and $M$ are given by symbolic
expressions with nested square roots of about 7000 operations in total,
and the task is to prove $N = -(1 - |M|^2)^2$. The program parses the
expressions with `gr_set_str`, computes $E = -(1 - |M|^2)^2$ and tests
$E = N$ with `gr_equal`.

By default it uses `qqbar_t` arithmetic. Every intermediate result is
reduced to its minimal polynomial, so the evaluation of $N$ dominates
(a few seconds), while the final comparison of two numbers of degree 8
is instantaneous. With the flag `-ca`, it uses `ca_t` arithmetic. The
evaluation is then much faster, since it consists of field operations
in a field with ten generators (nested square roots built on $\sqrt 2$),
but the result for $N$ is a large rational function of these
generators, and the final zero test, which has to evaluate it as an
algebraic number, takes most of the time. The total times of the two
approaches are similar, which is typical: the hard part of the problem
is moved, not removed.

The printed numerical value of $N$ in the `-ca` run is an indeterminate
enclosure, shown as `[+/- inf]` (as `nan` before FLINT 3.7): `ca_print`
evaluates the representation once at a fixed low precision, and the
denominator of the large rational function cannot be separated from
zero at that precision (see [](ca#printing)). The value itself is
correct, as the zero test shows.
