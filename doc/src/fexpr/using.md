---
title: Building and printing expressions
summary: How to construct symbolic expressions from atoms and calls, print them as text and LaTeX, inspect them with views, evaluate them numerically, replace subexpressions and expand arithmetic.
order: 1
---

This guide builds the expression for Euler's solution of the Basel problem, $\sum_{n=1}^\infty
1/n^2 = \pi^2/6$, step by step and then inspects, evaluates and
transforms it. It uses the builtin symbols `Sum`, `For`, `Infinity`,
`Pi`, `Equal`, `Div` and `Pow` of [](fexpr_builtin).

## Atoms and calls {#atoms}

An `fexpr_t` must be initialised with `fexpr_init` and cleared with
`fexpr_clear`, like other FLINT types. Atoms are created with
`fexpr_set_si`, `fexpr_set_fmpz` and `fexpr_set_fmpq` (integers and
fractions), `fexpr_set_symbol_str` (symbols by name, which gives a
builtin symbol if the name is builtin) or `fexpr_set_symbol_builtin`
(builtin symbols by their C constant, such as `FEXPR_Pi`), and
`fexpr_set_string`.

Calls are created with `fexpr_call0`, ..., `fexpr_call4` and
`fexpr_call_vec`, which take the function as an expression, with
`fexpr_call_builtin1` and `fexpr_call_builtin2`, which take a builtin
symbol constant, and with the shortcuts `fexpr_add`, `fexpr_mul`,
`fexpr_div`, `fexpr_pow` and so on for arithmetic. None of these
simplify anything: `fexpr_add` of `x` and `0` is `Add(x, 0)`. The
output of `fexpr_call0` to `fexpr_call_vec` must not be one of the
inputs; `fexpr_call_builtin1` and `fexpr_call_builtin2` allow it.

Since an expression is a flat array of words, every call copies its
arguments. Building a large expression by repeatedly wrapping it, for
example adding terms one at a time with binary `Add`, therefore takes
quadratic time; collect the terms in an array and call `fexpr_call_vec`
once instead.

## A first program {#program}

```c example name=fexpr_using
#include <flint/fexpr.h>
#include <flint/fexpr_builtin.h>

int main(void)
{
    fexpr_t n, one, two, six, inf, t, u, sum, rhs, eq, view;
    char * s;
    slong i;

    fexpr_init(n); fexpr_init(one); fexpr_init(two); fexpr_init(six);
    fexpr_init(inf); fexpr_init(t); fexpr_init(u); fexpr_init(sum);
    fexpr_init(rhs); fexpr_init(eq);

    fexpr_set_symbol_str(n, "n");
    fexpr_set_si(one, 1);
    fexpr_set_si(two, 2);
    fexpr_set_si(six, 6);
    fexpr_set_symbol_builtin(inf, FEXPR_Infinity);

    /* Sum(1/n^2, For(n, 1, Infinity)) */
    fexpr_pow(t, n, two);
    fexpr_div(u, one, t);
    fexpr_set_symbol_builtin(t, FEXPR_For);
    fexpr_call3(rhs, t, n, one, inf);
    fexpr_call_builtin2(sum, FEXPR_Sum, u, rhs);

    /* Pi^2 / 6 */
    fexpr_set_symbol_builtin(t, FEXPR_Pi);
    fexpr_pow(u, t, two);
    fexpr_div(rhs, u, six);

    fexpr_call_builtin2(eq, FEXPR_Equal, sum, rhs);

    fexpr_print(eq);
    flint_printf("\n");
    s = fexpr_get_str_latex(eq, 0);
    flint_printf("%s\n", s);
    flint_free(s);

    /* walk over the arguments of the Equal call with views */
    fexpr_view_arg(view, eq, 0);
    for (i = 0; i < fexpr_nargs(eq); i++)
    {
        flint_printf("argument %wd: depth %wd, %wd leaves\n", i,
            fexpr_depth(view), fexpr_num_leaves(view));
        fexpr_view_next(view);
    }

    /* numerical value of the right-hand side */
    s = fexpr_get_decimal_str(rhs, 20, 0);
    flint_printf("Pi^2/6 = %s\n", s);
    flint_free(s);

    /* replace Pi by x and expand (x + 1)^2 / 6 */
    fexpr_set_symbol_str(t, "x");
    fexpr_add(u, t, one);
    fexpr_set_symbol_builtin(t, FEXPR_Pi);
    fexpr_replace(rhs, rhs, t, u);
    fexpr_print(rhs);
    flint_printf("\n");
    fexpr_expanded_normal_form(rhs, rhs, 0);
    fexpr_print(rhs);
    flint_printf("\n");

    fexpr_clear(n); fexpr_clear(one); fexpr_clear(two); fexpr_clear(six);
    fexpr_clear(inf); fexpr_clear(t); fexpr_clear(u); fexpr_clear(sum);
    fexpr_clear(rhs); fexpr_clear(eq);
    return 0;
}
```
```output
Equal(Sum(Div(1, Pow(n, 2)), For(n, 1, Infinity)), Div(Pow(Pi, 2), 6))
\sum_{n=1}^{\infty} \frac{1}{{n}^{2}} = \frac{{\pi}^{2}}{6}
argument 0: depth 4, 10 leaves
argument 1: depth 3, 5 leaves
Pi^2/6 = 1.6449340668482264365
Div(Pow(Add(x, 1), 2), 6)
Div(Add(Pow(x, 2), Mul(2, x), 1), 6)
```

The LaTeX output renders as $\sum_{n=1}^{\infty} \frac{1}{{n}^{2}} =
\frac{{\pi}^{2}}{6}$.

## Inspecting expressions {#inspecting}

`fexpr_nargs` returns the number of arguments of a call (and $-1$ for an
atom), `fexpr_func` and `fexpr_arg` copy the function and the
arguments, and `fexpr_is_builtin_call` tests for a call of a given
builtin symbol, for example `fexpr_is_builtin_call(eq, FEXPR_Equal)`.

Copies are unnecessary for reading: `fexpr_view_func`,
`fexpr_view_arg` and `fexpr_view_next` set a *view*, an `fexpr_t` that
points into the data of the parent expression without owning it, as in
the loop above. A view is not initialised or cleared, must not be
modified, and becomes invalid when the parent is modified or cleared.

## Evaluation and transformation {#transforming}

Expressions have no evaluation semantics of their own. The module
offers a few transformations:

- `fexpr_get_decimal_str` evaluates an expression numerically in ball
  arithmetic, for the constants and functions it knows;
- `fexpr_replace` and its variants substitute subexpressions;
- `fexpr_expanded_normal_form` expands an expression as a rational
  function of its non-arithmetic subexpressions, which are treated as
  independent variables.

For exact evaluation, convert to `ca_t` with `ca_set_fexpr` (see
[](ca)) or to `qqbar_t` with `qqbar_set_fexpr` (see [](qqbar)); in the
other direction, `ca_get_fexpr`, `ca_poly_get_fexpr`,
`ca_mat_get_fexpr` and `qqbar_get_fexpr_formula` produce expressions,
which can then be printed as LaTeX.
