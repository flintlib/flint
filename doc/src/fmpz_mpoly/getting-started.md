---
title: Getting started with fmpz_mpoly
summary: How to set up a polynomial ring, choose a monomial ordering, build polynomials from strings, coefficients or terms, print them, iterate over their terms, and work with exponents of any size.
order: 1
---

An `fmpz_mpoly_t` is a sparse polynomial in several variables with
integer coefficients. Unlike the univariate types, it does not know its
ring: the number of variables and the monomial ordering are kept in a
context object, an `fmpz_mpoly_ctx_t`, which is passed to every function.
Everything here applies with small changes to `fmpq_mpoly`, `nmod_mpoly`,
`fmpz_mod_mpoly` and the other `mpoly` types.

## A first program {#first}

The program below works in $\Z[x, y, z]$. It parses two polynomials,
multiplies them, divides the product back, and looks up one coefficient.

```c example name=fmpz_mpoly_first
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>

int main(void)
{
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t f, g, h;
    const char * vars[] = {"x", "y", "z"};
    fmpz_t c;

    fmpz_mpoly_ctx_init(ctx, 3, ORD_LEX);    /* Z[x,y,z], x > y > z */
    fmpz_mpoly_init(f, ctx);
    fmpz_mpoly_init(g, ctx);
    fmpz_mpoly_init(h, ctx);
    fmpz_init(c);

    fmpz_mpoly_set_str_pretty(f, "x^2*y - 3*y*z + 2", vars, ctx);
    fmpz_mpoly_set_str_pretty(g, "(x + y + z)^2", vars, ctx);
    fmpz_mpoly_mul(h, f, g, ctx);

    flint_printf("f = "); fmpz_mpoly_print_pretty(f, vars, ctx); flint_printf("\n");
    flint_printf("g = "); fmpz_mpoly_print_pretty(g, vars, ctx); flint_printf("\n");
    flint_printf("f*g has %wd terms, total degree %wd, degree %wd in y\n",
        fmpz_mpoly_length(h, ctx), fmpz_mpoly_total_degree_si(h, ctx),
        fmpz_mpoly_degree_si(h, 1, ctx));

    fmpz_mpoly_divexact(h, h, g, ctx);       /* back to f */
    flint_printf("f*g/g == f: %d\n", fmpz_mpoly_equal(h, f, ctx));

    {
        ulong e[3] = {2, 1, 0};              /* the monomial x^2*y */
        fmpz_mpoly_get_coeff_fmpz_ui(c, f, e, ctx);
        flint_printf("coefficient of x^2*y in f: "); fmpz_print(c); flint_printf("\n");
    }

    fmpz_clear(c);
    fmpz_mpoly_clear(f, ctx);
    fmpz_mpoly_clear(g, ctx);
    fmpz_mpoly_clear(h, ctx);
    fmpz_mpoly_ctx_clear(ctx);
    return 0;
}
```
```output
f = x^2*y-3*y*z+2
g = x^2+2*x*y+2*x*z+y^2+2*y*z+z^2
f*g has 18 terms, total degree 5, degree 3 in y
f*g/g == f: 1
coefficient of x^2*y in f: 1
```

Polynomials must be cleared before their context. The context only
stores the number of variables and the ordering; the names of the
variables are not part of it and are passed to the input and output
functions when needed.

## Variables and monomial orderings {#orderings}

The variables are numbered from 0, and the variable of index 0 is the
most significant one in the ordering: with the names `{"x", "y", "z"}`
above, $x > y > z$. Exponent vectors such as `e` above are indexed in the
same way, and so are the arguments `var` of functions such as
`fmpz_mpoly_degree_si`.

The terms of a polynomial are stored in decreasing order for the monomial
ordering of the context, and printed in that order. Three orderings are
available:

- `ORD_LEX`: compare the exponents of $x$, then of $y$, and so on.
- `ORD_DEGLEX`: compare total degrees first, then lexicographically.
- `ORD_DEGREVLEX`: compare total degrees first; among monomials of the
  same degree, the one with the *smaller* exponent of the *last* variable
  in which they differ is larger.

```c example name=fmpz_mpoly_orderings
#include <flint/flint.h>
#include <flint/fmpz_mpoly.h>

int main(void)
{
    const char * vars[] = {"x", "y", "z"};
    const ordering_t ords[3] = {ORD_LEX, ORD_DEGLEX, ORD_DEGREVLEX};
    const char * names[3] = {"lex", "deglex", "degrevlex"};
    int k;

    for (k = 0; k < 3; k++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f;

        fmpz_mpoly_ctx_init(ctx, 3, ords[k]);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_set_str_pretty(f, "x*z^2 + y^3 + x^2 + x*y*z + z^3", vars, ctx);
        flint_printf("%-10s ", names[k]);
        fmpz_mpoly_print_pretty(f, vars, ctx);
        flint_printf("\n");
        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }
    return 0;
}
```
```output
lex        x^2+x*y*z+x*z^2+y^3+z^3
deglex     x*y*z+x*z^2+y^3+z^3+x^2
degrevlex  y^3+x*y*z+x*z^2+z^3+x^2
```

The ordering matters for division with remainder, for Gröbner bases and
for what "leading term" means, and it influences which multiplication
algorithm is fastest; for most other purposes any ordering will do. A
polynomial cannot change context; to move a polynomial to a ring with
other variables or another ordering, use
`fmpz_mpoly_compose_fmpz_mpoly_gen`.

## Parsing and printing {#strings}

`fmpz_mpoly_set_str_pretty` parses an expression built from nonnegative
integers, the variable names, `+`, `-`, `*`, `/` (exact division only),
`^` followed by a nonnegative integer, and parentheses. It returns 0 on
success and $-1$ if the string cannot be parsed. Implicit multiplication
(`2x`), negative exponents, unknown names and inexact divisions are
errors. If the names are `NULL`, the variables are called `x1`, `x2`, and
so on, both for parsing and for printing. `fmpz_mpoly_get_str_pretty`
returns a string that must be freed with `flint_free`.

```c example name=fmpz_mpoly_strings
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>

int main(void)
{
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t f;
    const char * vars[] = {"x", "y"};
    char * s;
    int ret;

    fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);
    fmpz_mpoly_init(f, ctx);

    ret = fmpz_mpoly_set_str_pretty(f, "(x - 2*y)^3 / (x - 2*y) + 1", vars, ctx);
    s = fmpz_mpoly_get_str_pretty(f, vars, ctx);
    flint_printf("%d: %s\n", ret, s);
    flint_free(s);

    s = fmpz_mpoly_get_str_pretty(f, NULL, ctx);    /* default names */
    flint_printf("default names: %s\n", s);
    flint_free(s);

    ret = fmpz_mpoly_set_str_pretty(f, "x1^2 - 3*x2", NULL, ctx);
    flint_printf("%d: ", ret);
    fmpz_mpoly_print_pretty(f, vars, ctx);
    flint_printf("\n");

    flint_printf("2x: %d\n", fmpz_mpoly_set_str_pretty(f, "2x", vars, ctx));
    flint_printf("x/y: %d\n", fmpz_mpoly_set_str_pretty(f, "x/y", vars, ctx));
    flint_printf("x^-1: %d\n", fmpz_mpoly_set_str_pretty(f, "x^-1", vars, ctx));
    flint_printf("z: %d\n", fmpz_mpoly_set_str_pretty(f, "z", vars, ctx));

    fmpz_mpoly_clear(f, ctx);
    fmpz_mpoly_ctx_clear(ctx);
    return 0;
}
```
```output
0: x^2-4*x*y+4*y^2+1
default names: x1^2-4*x1*x2+4*x2^2+1
0: x^2-3*y
2x: -1
x/y: -1
x^-1: -1
z: -1
```

After a failed parse the polynomial holds an unspecified (but valid)
value. Parsing is convenient but slow for large inputs: it evaluates the
expression with polynomial arithmetic. To build large polynomials, use
the term functions of the next section.

## Building a polynomial term by term {#building}

There are three ways to set coefficients:

- `fmpz_mpoly_set_coeff_fmpz_ui` and its variants set the coefficient of
  a given monomial, keeping the polynomial sorted. Each call costs a
  binary search plus, when a term is inserted or removed, a move of the
  following terms, so building a polynomial of $n$ terms this way can
  take $O(n^2)$ time.
- `fmpz_mpoly_push_term_fmpz_ui` and its variants append a term at the
  end in amortised constant time without any checks. If the terms are
  appended in decreasing order with distinct monomials and nonzero
  coefficients, the result is directly valid.
- If the terms come in no particular order, push them all and then call
  `fmpz_mpoly_sort_terms` and `fmpz_mpoly_combine_like_terms`, which
  merges repeated monomials and removes zero coefficients. This is the
  fast way to build a large polynomial.

The following program builds $\sum_{0 \le i, j < 3} (i - j) x^i y^j +
5x^2y^2$ in degrevlex order from unsorted terms, and then reads it back
term by term.

```c example name=fmpz_mpoly_terms
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>

int main(void)
{
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t f;
    const char * vars[] = {"x", "y"};
    ulong e[2];
    slong i, j;
    fmpz_t c;

    fmpz_mpoly_ctx_init(ctx, 2, ORD_DEGREVLEX);
    fmpz_mpoly_init(f, ctx);
    fmpz_init(c);

    /* sum of (i - j) x^i y^j for 0 <= i, j < 3, terms appended in any order */
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
        {
            e[0] = i;
            e[1] = j;
            fmpz_mpoly_push_term_si_ui(f, i - j, e, ctx);
        }
    fmpz_mpoly_push_term_si_ui(f, 5, e, ctx);    /* another x^2 y^2 */

    flint_printf("canonical before: %d\n", fmpz_mpoly_is_canonical(f, ctx));
    fmpz_mpoly_sort_terms(f, ctx);
    fmpz_mpoly_combine_like_terms(f, ctx);       /* also drops zero terms */
    flint_printf("canonical after:  %d\n", fmpz_mpoly_is_canonical(f, ctx));
    flint_printf("f = "); fmpz_mpoly_print_pretty(f, vars, ctx); flint_printf("\n");

    for (i = 0; i < fmpz_mpoly_length(f, ctx); i++)
    {
        fmpz_mpoly_get_term_coeff_fmpz(c, f, i, ctx);
        fmpz_mpoly_get_term_exp_ui(e, f, i, ctx);
        flint_printf("term %wd: coefficient ", i);
        fmpz_print(c);
        flint_printf(", exponents (%wu, %wu)\n", e[0], e[1]);
    }

    fmpz_clear(c);
    fmpz_mpoly_clear(f, ctx);
    fmpz_mpoly_ctx_clear(ctx);
    return 0;
}
```
```output
canonical before: 0
canonical after:  1
f = 5*x^2*y^2+x^2*y-x*y^2+2*x^2-2*y^2+x-y
term 0: coefficient 5, exponents (2, 2)
term 1: coefficient 1, exponents (2, 1)
term 2: coefficient -1, exponents (1, 2)
term 3: coefficient 2, exponents (2, 0)
term 4: coefficient -2, exponents (0, 2)
term 5: coefficient 1, exponents (1, 0)
term 6: coefficient -1, exponents (0, 1)
```

The terms $(i - i) x^i y^i$ for $i = 0, 1$ have zero coefficient and
disappear, and the two $x^2 y^2$ terms are merged. Until the polynomial
is canonical again, it must not be passed to any function except the
container functions (see [](fmpz_mpoly#conventions)): arithmetic
functions assume sorted terms and nonzero coefficients and give wrong
results otherwise.

## Iterating over terms {#iterating}

The terms of a polynomial are indexed from 0 (the leading term) to
`fmpz_mpoly_length(f, ctx) - 1`. For each term there are functions to
read

- the coefficient: `fmpz_mpoly_get_term_coeff_fmpz`, or the `_ui` and
  `_si` versions when it is known to fit, or a pointer to it with
  `fmpz_mpoly_term_coeff_ref`;
- the exponent vector: `fmpz_mpoly_get_term_exp_ui`, `_si` or `_fmpz`;
- a single exponent: `fmpz_mpoly_get_term_var_exp_ui`;
- the term or its monomial as a polynomial: `fmpz_mpoly_get_term` and
  `fmpz_mpoly_get_term_monomial`.

Each access costs time proportional to the number of variables. To see a polynomial as a
polynomial in one variable with multivariate coefficients, convert it
with `fmpz_mpoly_to_univar`; to read the coefficient of a power product
of some of the variables, use `fmpz_mpoly_get_coeff_vars_ui`.

## Exponents of any size {#exponents}

Exponents are not limited to a machine word. Internally, the exponents of
a monomial are packed into bit fields of a common width, the field `bits`
of the polynomial, which every function enlarges as needed; the details
are described in [](algorithms/multivariate/representation). Exponent
vectors of type `fmpz` give access to exponents of any size, and the
`_ui` and `_si` access functions throw an exception if an exponent does
not fit.

The program below shows the packing width growing as larger exponents
appear, and the behaviour of powering with a huge exponent:
`fmpz_mpoly_pow_fmpz` returns 0 instead of attempting a computation whose
result cannot be stored, which is possible only for monomials with
coefficient $\pm 1$.

```c example name=fmpz_mpoly_bigexp
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mpoly.h>

int main(void)
{
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t f, g;
    const char * vars[] = {"x", "y"};
    fmpz_t e0, e1, k;
    fmpz * exp[2];

    fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);
    fmpz_mpoly_init(f, ctx);
    fmpz_mpoly_init(g, ctx);
    fmpz_init(e0); fmpz_init(e1); fmpz_init(k);
    exp[0] = e0;
    exp[1] = e1;

    fmpz_mpoly_set_str_pretty(f, "x^1000*y + 1", vars, ctx);
    flint_printf("bits = %wu\n", f->bits);

    fmpz_one(e0);
    fmpz_mul_2exp(e0, e0, 62);                  /* x^(2^62) */
    fmpz_zero(e1);
    fmpz_mpoly_set_coeff_si_fmpz(f, 3, exp, ctx);
    flint_printf("bits = %wu\n", f->bits);

    fmpz_mpoly_mul(g, f, f, ctx);               /* contains x^(2^63) */
    flint_printf("bits = %wu, terms = %wd\n", g->bits, fmpz_mpoly_length(g, ctx));
    flint_printf("exponents of term 0 fit ulong: %d, slong: %d\n",
                 fmpz_mpoly_term_exp_fits_ui(g, 0, ctx),
                 fmpz_mpoly_term_exp_fits_si(g, 0, ctx));
    fmpz_mpoly_get_term_exp_fmpz(exp, g, 0, ctx);
    flint_printf("leading monomial x^"); fmpz_print(e0);
    flint_printf(" y^"); fmpz_print(e1); flint_printf("\n");

    fmpz_set_str(k, "1000000000000000000000", 10);   /* 10^21 */
    fmpz_mpoly_set_str_pretty(f, "x + 1", vars, ctx);
    flint_printf("(x + 1)^k: %d\n", fmpz_mpoly_pow_fmpz(g, f, k, ctx));
    fmpz_mpoly_set_str_pretty(f, "-x*y^2", vars, ctx);
    flint_printf("(-x*y^2)^k: %d, ", fmpz_mpoly_pow_fmpz(g, f, k, ctx));
    fmpz_mpoly_print_pretty(g, vars, ctx);
    flint_printf("\n");

    fmpz_clear(e0); fmpz_clear(e1); fmpz_clear(k);
    fmpz_mpoly_clear(f, ctx);
    fmpz_mpoly_clear(g, ctx);
    fmpz_mpoly_ctx_clear(ctx);
    return 0;
}
```
```output
bits = 32
bits = 64
bits = 128, terms = 6
exponents of term 0 fit ulong: 1, slong: 0
leading monomial x^9223372036854775808 y^0
(x + 1)^k: 0
(-x*y^2)^k: 1, x^1000000000000000000000*y^2000000000000000000000
```

On a 64-bit machine, the two exponents of $x^{1000} y$ are packed into one
word with 32 bits each, the largest width that still fits both fields in
one word. The exponent $2^{62}$ needs 64 bits per field (one bit is kept
free to detect overflow), and $2^{63}$ needs two words per field. Wider
packings make every monomial operation more expensive, so performance is
best when the exponents stay small.

## Conventions worth knowing {#conventions}

- **Aliasing.** Outputs may be the same objects as inputs, as in
  `fmpz_mpoly_divexact(h, h, g, ctx)` above. Several outputs of one
  function (quotient and remainder, GCD and cofactors) must be distinct.
- **Failure.** Functions that can fail return an `int`, 1 for success.
  Besides a few algorithm-specific functions, failure means that the
  result would be unreasonably large, as in the powering example; the
  general purpose functions such as `fmpz_mpoly_gcd` do not fail in
  practice, but their return value should still be checked.
- **Exactness.** `fmpz_mpoly_divides` tests divisibility and returns the
  quotient if it exists; `fmpz_mpoly_divexact` throws if the division is
  not exact; `fmpz_mpoly_div` and `fmpz_mpoly_divrem` perform division
  with remainder, which over $\Z$ is subtle, see [](fmpz_mpoly/division).
- **Threads.** Multiplication and exact division use several threads for
  large inputs if FLINT is allowed to, see `flint_set_num_threads`.

The pages [](fmpz_mpoly/mul), [](fmpz_mpoly/division) and
[](fmpz_mpoly/gcd) describe the main operations and their algorithms,
factorisation is in [](fmpz_mpoly_factor), rational functions in
[](fmpz_mpoly_q), and the algorithms of the whole family of `mpoly`
modules are compared in [](algorithms/multivariate/overview).
