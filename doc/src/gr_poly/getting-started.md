---
title: Getting started with gr_poly
summary: How to create and use polynomials over a generic ring, check status codes, write code that works over any ring, and deal with coefficients whose zero test is undecidable.
order: 1
---

A `gr_poly_t` is a dense polynomial whose coefficients are elements of a
ring described by a context `gr_ctx_t` ([](gr)). The same functions work
over the integers, the rationals, $\Z/n\Z$, finite fields, real and
complex balls, floating-point numbers and any other ring with a
context.

## A first program {#first}

Every function takes the context of the coefficient ring as its last
argument. Polynomials are initialised with `gr_poly_init`, can be set
from a string with `gr_poly_set_str`, and are cleared with
`gr_poly_clear`; contexts are cleared with `gr_ctx_clear`. The program
below divides $x^4 - 1$ by $2x^2 + 3x + 1$ over $\mathbb{Q}$ and computes a GCD.

```c example name=gr_poly_first
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_poly.h>

int main(void)
{
    gr_ctx_t QQ;
    gr_poly_t f, g, q, r, d;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpq(QQ);
    gr_poly_init(f, QQ);
    gr_poly_init(g, QQ);
    gr_poly_init(q, QQ);
    gr_poly_init(r, QQ);
    gr_poly_init(d, QQ);

    status |= gr_poly_set_str(f, "x^4 - 1", "x", QQ);
    status |= gr_poly_set_str(g, "2*x^2 + 3*x + 1", "x", QQ);

    status |= gr_poly_divrem(q, r, f, g, QQ);
    status |= gr_poly_gcd(d, f, g, QQ);

    flint_printf("q = "); status |= gr_poly_print(q, QQ); flint_printf("\n");
    flint_printf("r = "); status |= gr_poly_print(r, QQ); flint_printf("\n");
    flint_printf("gcd = "); status |= gr_poly_print(d, QQ); flint_printf("\n");
    flint_printf("status = %d\n", status);

    gr_poly_clear(f, QQ);
    gr_poly_clear(g, QQ);
    gr_poly_clear(q, QQ);
    gr_poly_clear(r, QQ);
    gr_poly_clear(d, QQ);
    gr_ctx_clear(QQ);
    return 0;
}
```
```output
q = (7/8) + (-3/4)*x + (1/2)*x^2
r = (-15/8) + (-15/8)*x
gcd = 1 + x
status = 0
```

Polynomials are printed in increasing order of degree. The GCD over a
field is monic. The status codes of all calls are combined with `|`,
the usual idiom: `GR_SUCCESS` is 0, so the combined value is zero
exactly when every call succeeded. Most functions are declared with
`WARN_UNUSED_RESULT`, so the compiler reminds you to look at them; in
test code and examples, `GR_MUST_SUCCEED(...)` aborts the program on
failure.

## Code for any ring {#generic}

A function that takes a context works over every ring. The function
`show_gcd` below computes $\gcd(2x^3 + 6x^2 - 2x - 6, 4x^2 + 8x - 12)$,
and `main` calls it over $\Z$, $\mathbb{Q}$, $\Z/5\Z$ and $\Z/6\Z$.

```c example name=gr_poly_generic_gcd
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_poly.h>

/* The same code for any coefficient ring. */
static void
show_gcd(gr_ctx_t ctx)
{
    gr_poly_t f, g, d;
    int status = GR_SUCCESS;

    gr_poly_init(f, ctx);
    gr_poly_init(g, ctx);
    gr_poly_init(d, ctx);

    status |= gr_poly_set_str(f, "2*x^3 + 6*x^2 - 2*x - 6", "x", ctx);
    status |= gr_poly_set_str(g, "4*x^2 + 8*x - 12", "x", ctx);
    status |= gr_poly_gcd(d, f, g, ctx);

    status |= gr_ctx_print(ctx);
    flint_printf(": ");
    if (status == GR_SUCCESS)
        status |= gr_poly_print(d, ctx);
    else
        flint_printf("failed (status %d)", status);
    flint_printf("\n");

    gr_poly_clear(f, ctx);
    gr_poly_clear(g, ctx);
    gr_poly_clear(d, ctx);
}

int main(void)
{
    gr_ctx_t ctx;

    gr_ctx_init_fmpz(ctx);
    show_gcd(ctx);
    gr_ctx_clear(ctx);

    gr_ctx_init_fmpq(ctx);
    show_gcd(ctx);
    gr_ctx_clear(ctx);

    gr_ctx_init_nmod(ctx, 5);
    show_gcd(ctx);
    gr_ctx_clear(ctx);

    gr_ctx_init_nmod(ctx, 6);
    show_gcd(ctx);
    gr_ctx_clear(ctx);
    return 0;
}
```
```output
Integer ring (fmpz): -6 + 4*x + 2*x^2
Rational field (fmpq): -3 + 2*x + x^2
Integers mod 5 (_gr_nmod): 2 + 2*x + x^2
Integers mod 6 (_gr_nmod): failed (status 2)
```

Over $\Z$ the GCD keeps the common
content 2 and is normalised to a positive leading coefficient; over the
fields $\mathbb{Q}$ and $\Z/5\Z$ it is monic. $(\Z/6\Z)[x]$ is not a GCD domain,
and the function reports that it cannot compute an answer instead of
returning a wrong one. Behind the same call, the algorithms differ too:
the subresultant algorithm over $\Z$, the Euclidean algorithm over $\mathbb{Q}$
and the GCD of `nmod_poly` over $\Z/5\Z$ ([](gr_poly/gcd)).

## Status codes {#status}

Every function that can fail returns a status, a combination of two
flags:

- `GR_DOMAIN` means that the result does not exist in the ring:
  division by zero, a division that is not possible in the ring (the
  quotient of $x^2$ by $2x$ in $\Z[x]$), the power series inverse of a
  series whose constant term is not invertible, or the square root of a
  series whose constant term is not a square.
- `GR_UNABLE` means that the result may exist but could not be
  computed: a zero test that cannot be decided, a case that is not
  implemented, a ring that does not have the required structure (as
  for the GCD over $\Z/6\Z$ above), or a result too large to represent.

A function that fails leaves its outputs in an undefined but valid
state: they can be cleared or overwritten, but their values are
meaningless. Functions that answer a yes/no question return a
`truth_t`: `T_TRUE`, `T_FALSE` or `T_UNKNOWN`.

## Inexact coefficients and weak normalisation {#inexact}

Over rings such as real balls (`gr_ctx_init_real_arb`), a coefficient
can be an enclosure like $[0 \pm 10^{-10}]$ that may or may not be zero.
`gr_poly` normalises polynomials by removing leading coefficients only
when `gr_is_zero` returns `T_TRUE`, so such a coefficient stays: the
length of the polynomial is then only an upper bound for its degree
plus one. Functions whose result depends on the exact degree check the
leading coefficient and return `GR_UNABLE` when it is not provably
nonzero:

```c example name=gr_poly_inexact
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_poly.h>

int main(void)
{
    gr_ctx_t RR;
    gr_poly_t f, g, h, q, r;

    gr_ctx_init_real_arb(RR, 53);
    gr_poly_init(f, RR);
    gr_poly_init(g, RR);
    gr_poly_init(h, RR);
    gr_poly_init(q, RR);
    gr_poly_init(r, RR);

    GR_MUST_SUCCEED(gr_poly_set_str(f, "x^5 + 1", "x", RR));
    /* g = x^2 + 3x + 1 with an uncertain leading coefficient */
    GR_MUST_SUCCEED(gr_poly_set_str(g, "[1 +/- 1e-10]*x^2 + 3*x + 1", "x", RR));
    flint_printf("divrem by g:  status %d\n", gr_poly_divrem(q, r, f, g, RR));

    /* h = g - x^2 has the leading coefficient [+/- 1e-10] */
    GR_MUST_SUCCEED(gr_poly_set_str(h, "x^2", "x", RR));
    GR_MUST_SUCCEED(gr_poly_sub(h, g, h, RR));
    flint_printf("h = "); gr_poly_print(h, RR); flint_printf("\n");
    flint_printf("length of h:  %wd\n", gr_poly_length(h, RR));
    flint_printf("divrem by h:  status %d\n", gr_poly_divrem(q, r, f, h, RR));

    gr_poly_clear(f, RR);
    gr_poly_clear(g, RR);
    gr_poly_clear(h, RR);
    gr_poly_clear(q, RR);
    gr_poly_clear(r, RR);
    gr_ctx_clear(RR);
    return 0;
}
```
```output
divrem by g:  status 0
h = 1 + 3.000000000000000*x + [+/- 1.01e-10]*x^2
length of h:  3
divrem by h:  status 2
```

Dividing by $g$ works because its leading coefficient $[1 \pm 10^{-10}]$
excludes zero. The leading coefficient of $h$ contains zero, so the
quotient by $h$ could have degree 3 or 4 (or $h$ could have degree 1),
and `gr_poly_divrem` returns `GR_UNABLE` (2). Operations that do not
depend on the degree, such as addition, multiplication, evaluation and
power series arithmetic, work as usual and produce enclosures of the
exact results. The same rules hold for the exact real and complex
numbers of `ca`, where a zero test may be undecidable, and the design is
described in [](algorithms/generic/overview).

## Low-level functions {#underscore}

Functions with a leading underscore work on coefficient arrays: a
pointer of type `gr_ptr` (or `gr_srcptr` for inputs) and an explicit
length. They do not allocate the output, do not normalise it, and
usually require nonzero lengths and distinct input and output arrays.
The element at index $i$ of an array `v` is `GR_ENTRY(v, i, sz)` with
`sz = ctx->sizeof_elem`, and the coefficients of a `gr_poly_t` are
`poly->coeffs` with length `poly->length`. To call an underscore
function on polynomials, allocate the output with `gr_poly_fit_length`
and finish with `_gr_poly_set_length_normalise`, as the non-underscore
functions do. Temporary arrays of elements are created with
`GR_TMP_INIT_VEC` and freed with `GR_TMP_CLEAR_VEC`.

## Polynomials as ring elements {#poly-ring}

The context of a polynomial ring $R[x]$, created by
`gr_ctx_init_gr_poly(Rx, R)`, has `gr_poly_struct` elements. It is
needed to build structures over $R[x]$: matrices of polynomials
(`gr_mat`), polynomials over polynomials (bivariate polynomials), or
power series over polynomials. For arithmetic on single polynomials, the functions of
this module with the context of $R$ are simpler and faster.
