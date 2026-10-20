---
title: Getting started with generic rings
summary: Contexts, elements, status codes and truth values, vectors, printing and conversions, with programs that run the same code over several rings.
order: 1
---

This guide introduces the `gr` interface through small complete
programs: one generic function evaluated over six rings, the ways of
allocating elements, the meaning of status codes and truth values,
vectors, printing and parsing, and conversions between rings. All
functions mentioned here are documented on the module page [`gr`](gr)
and in [`gr_vec`](gr_vec); the constructors of the rings are listed in
[`gr_domains`](gr_domains).

## A first program {#first}

A generic function takes its elements as `gr_ptr` (outputs) and
`gr_srcptr` (inputs) and the ring as a final `gr_ctx_t` argument, and it
returns a status code. The function below computes the partial sum
$\sum_{k=0}^{n} x^k / k!$ of the exponential series. It is evaluated at
$x = 1$ in $\Z$, $\mathbb{Q}$, $\mathbb{R}$ (balls with 64-bit precision),
$\Z/7\Z$ and $\Z/10\Z$, and at the generator $x$ of $\mathbb{Q}[x]$.

```c example name=gr_first
#include <flint/flint.h>
#include <flint/gr.h>

/* res = sum_{k=0}^{n} x^k / k!, computed in any ring */
int exp_partial_sum(gr_ptr res, gr_srcptr x, ulong n, gr_ctx_t ctx)
{
    gr_ptr t;
    ulong k;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);
    status |= gr_one(t, ctx);            /* t = x^k / k! */
    status |= gr_one(res, ctx);
    for (k = 1; k <= n; k++)
    {
        status |= gr_mul(t, t, x, ctx);
        status |= gr_div_ui(t, t, k, ctx);
        status |= gr_add(res, res, t, ctx);
    }
    GR_TMP_CLEAR(t, ctx);
    return status;
}

void run(gr_ctx_t ctx, const char * xstr)
{
    gr_ptr x, y;
    int status;

    x = gr_heap_init(ctx);
    y = gr_heap_init(ctx);
    status = gr_set_str(x, xstr, ctx);
    status |= exp_partial_sum(y, x, 6, ctx);
    gr_ctx_print(ctx);
    if (status == GR_SUCCESS)
        flint_printf(": %{gr}\n", y, ctx);
    else
        flint_printf(": status %d\n", status);
    gr_heap_clear(x, ctx);
    gr_heap_clear(y, ctx);
}

int main(void)
{
    gr_ctx_t ZZ, QQ, RR, ZZ7, ZZ10, QQx;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_fmpq(QQ);
    gr_ctx_init_real_arb(RR, 64);
    gr_ctx_init_nmod(ZZ7, 7);
    gr_ctx_init_nmod(ZZ10, 10);
    gr_ctx_init_gr_poly(QQx, QQ);

    run(ZZ, "1");
    run(QQ, "1");
    run(RR, "1");
    run(ZZ7, "1");
    run(ZZ10, "1");
    run(QQx, "x");

    gr_ctx_clear(QQx);
    gr_ctx_clear(ZZ);
    gr_ctx_clear(QQ);
    gr_ctx_clear(RR);
    gr_ctx_clear(ZZ7);
    gr_ctx_clear(ZZ10);
    return 0;
}
```
```output
Integer ring (fmpz): status 1
Rational field (fmpq): 1957/720
Real numbers (arb, prec = 64): [2.718055555555555555 +/- 8.99e-19]
Integers mod 7 (_gr_nmod): 3
Integers mod 10 (_gr_nmod): status 1
Ring of polynomials over Rational field (fmpq): 1 + x + (1/2)*x^2 + (1/6)*x^3 + (1/24)*x^4 + (1/120)*x^5 + (1/720)*x^6
```

In $\Z$ the first division $1/2$ has no result, so the status is
`GR_DOMAIN` (1); in $\Z/10\Z$ the same happens because $2$ is not
invertible modulo 10. Modulo 7 every $k \le 6$ is invertible and
$1957/720 \equiv 4 / 6 \equiv 3$. Over the reals the result is a ball
that contains the exact rational value. The polynomial ring $\mathbb{Q}[x]$ is
built on top of the context `QQ`, and must be cleared before it.

## Contexts {#contexts}

Every ring is a `gr_ctx_t`, initialised by a constructor and freed with
`gr_ctx_clear`. Some common rings:

| Ring | Constructor | Element type |
|---|---|---|
| $\Z$, $\mathbb{Q}$, $\Z[i]$ | `gr_ctx_init_fmpz`, `gr_ctx_init_fmpq`, `gr_ctx_init_fmpzi` | `fmpz`, `fmpq`, `fmpzi_struct` |
| $\Z/n\Z$ | `gr_ctx_init_nmod` ($n < 2^{\bits}$), `gr_ctx_init_mpn_mod` ($n < 2^{1024}$), `gr_ctx_init_fmpz_mod` (any $n$) | `ulong`, limb array, `fmpz` |
| $\mathbb{F}_{p^d}$ | `gr_ctx_init_fq`, `gr_ctx_init_fq_nmod`, `gr_ctx_init_fq_zech` | `fq_struct`, ... |
| $\mathbb{R}$, $\mathbb{C}$ (balls) | `gr_ctx_init_real_arb`, `gr_ctx_init_complex_acb` | `arb_struct`, `acb_struct` |
| $\overline{\mathbb{Q}}$ | `gr_ctx_init_real_qqbar`, `gr_ctx_init_complex_qqbar` | `qqbar_struct` |
| exact reals and complex numbers | `gr_ctx_init_real_ca`, `gr_ctx_init_complex_ca` | `ca_struct` |
| $R[x]$, $R[x_1,\ldots,x_n]$, $R[[x]]$ | `gr_ctx_init_gr_poly`, `gr_ctx_init_gr_mpoly`, `gr_series_ctx_init` | `gr_poly_struct`, ... |
| $M_n(R)$, $\operatorname{Frac}(R)$, $R[x]/(m)$ | `gr_ctx_init_matrix_ring`, `gr_ctx_init_gr_fraction`, `gr_ctx_init_gr_poly_quotient` | `gr_mat_struct`, ... |

The full list, including floating-point contexts, number fields and
groups, is on [`gr_domains`](gr_domains). A context is cheap to create:
most constructors allocate nothing and only fill in a few words and a
pointer to a method table that is built once per type, so that
algorithms can create thousands of contexts, one per prime of a
multimodular computation ([](algorithms/generic/method-tables#cheap)).

A context answers questions about its ring: `gr_ctx_print` and
`%{gr_ctx}` describe it, and predicates such as `gr_ctx_is_field` and
`gr_ctx_is_exact` return truth values that generic algorithms use to
choose a method. For example `gr_ctx_is_field` is `T_TRUE` for
`gr_ctx_init_nmod(ctx, 7)`, `T_FALSE` for modulus 6, and `T_UNKNOWN` for
`gr_ctx_init_fmpz_mod` with a large modulus until primality is declared
with `gr_ctx_set_is_field`.

## Elements {#elements}

An element is `ctx->sizeof_elem` bytes of memory holding a value of the
ring's element type; it must be initialised (to zero) before use and
cleared afterwards. There are three ways to allocate one:

* `GR_TMP_INIT(x, ctx)` ... `GR_TMP_CLEAR(x, ctx)` (and the variants for
  two to five elements) allocate on the stack. This is the fastest way
  and the usual one inside functions, but the space is only released
  when the function returns, so never use these macros in a loop.
* `x = gr_heap_init(ctx)` ... `gr_heap_clear(x, ctx)` allocate on the
  heap; `gr_heap_init_vec` and `gr_heap_clear_vec` do the same for
  arrays.
* When the ring is known, a variable of the element type can be used
  directly: an `fmpz_t` is a valid element of `gr_ctx_init_fmpz`, an
  `arb_t` of `gr_ctx_init_real_arb`, a `ulong` of `gr_ctx_init_nmod`.
  This is how type-specific code calls generic algorithms.

Elements carry no reference to their ring. The caller must always pass
the context that the element belongs to; passing the wrong one is
undefined behaviour, just as with `fmpz_mod` and its context. Values are
never converted implicitly: to use an integer in a polynomial ring, call
a conversion such as `gr_set_other`, `gr_set_fmpz` or a mixed function
such as `gr_mul_fmpz`.

```c example name=gr_elements
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/arb.h>
#include <flint/gr.h>
#include <flint/gr_special.h>

int main(void)
{
    gr_ctx_t ZZ, RR, F7, F7x, ZZx;
    fmpz_t a;
    arb_t b;
    gr_ptr f, g;
    int status = GR_SUCCESS;

    /* elements of a known type can be declared with that type */
    gr_ctx_init_fmpz(ZZ);
    fmpz_init_set_ui(a, 3);
    status |= gr_pow_ui(a, a, 100, ZZ);
    flint_printf("3^100 = %{fmpz}\n", a);

    gr_ctx_init_real_arb(RR, 64);
    arb_init(b);
    status |= gr_set_fmpz(b, a, RR);
    status |= gr_log(b, b, RR);
    flint_printf("log(3^100) = %{gr}\n", b, RR);

    /* polynomials over Z/7Z, and a conversion from Z[x] */
    gr_ctx_init_nmod(F7, 7);
    gr_ctx_init_gr_poly(F7x, F7);
    gr_ctx_init_fmpz_poly(ZZx);
    f = gr_heap_init(F7x);
    g = gr_heap_init(ZZx);
    status |= gr_set_str(f, "(x + 1)^7", F7x);
    flint_printf("(x + 1)^7 = %{gr} in %{gr_ctx}\n", f, F7x, F7x);
    status |= gr_set_str(g, "10*x^2 + 15*x - 1", ZZx);
    status |= gr_set_other(f, g, ZZx, F7x);
    flint_printf("10*x^2 + 15*x - 1 = %{gr}\n", f, F7x);
    flint_printf("element sizes: %wd %wd %wd %wd\n", ZZ->sizeof_elem,
        RR->sizeof_elem, F7->sizeof_elem, F7x->sizeof_elem);

    gr_heap_clear(f, F7x);
    gr_heap_clear(g, ZZx);
    gr_ctx_clear(ZZx);
    gr_ctx_clear(F7x);
    gr_ctx_clear(F7);
    arb_clear(b);
    fmpz_clear(a);
    gr_ctx_clear(RR);
    gr_ctx_clear(ZZ);
    flint_printf("status %d\n", status);
    return 0;
}
```
```output
3^100 = 515377520732011331036461129765621272702107522001
log(3^100) = [109.8612288668109691 +/- 4.18e-17]
(x + 1)^7 = 1 + x^7 in Ring of polynomials over Integers mod 7 (_gr_nmod)
10*x^2 + 15*x - 1 = 6 + x + 3*x^2
element sizes: 8 48 8 24
status 0
```

The element of `F7x` is a `gr_poly_struct` (24 bytes: a pointer to the
coefficients, the allocation and the length), whose coefficients are
`ulong` residues. In a context built on a base ring, the element of the
outer ring owns memory allocated through the base ring, so the base ring
must outlive all elements of the outer ring.

## Status codes and truth values {#status}

Every function that can fail returns a status: `GR_SUCCESS` (0),
`GR_DOMAIN` (1) when the result is not defined in the ring, or
`GR_UNABLE` (2) when the implementation cannot compute it. The two
flags are combined with `|`, and checking once after a sequence of
operations is the usual pattern, as in `exp_partial_sum` above. When
the status is not `GR_SUCCESS`, the outputs hold unspecified (but valid)
values. `GR_MUST_SUCCEED(expr)` aborts on failure, for code in which
failure is impossible.

Predicates such as `gr_equal`, `gr_is_zero` and `gr_ctx_is_field`
return a `truth_t`: `T_TRUE`, `T_FALSE` or `T_UNKNOWN`. Since `T_TRUE`
is 0, never use a `truth_t` as a C condition; compare it with `T_TRUE`
or `T_FALSE`. The program below shows typical outcomes.

```c example name=gr_status
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_special.h>

const char * status_str(int status)
{
    if (status == GR_SUCCESS) return "GR_SUCCESS";
    if (status == GR_DOMAIN) return "GR_DOMAIN";
    if (status == GR_UNABLE) return "GR_UNABLE";
    return "GR_DOMAIN | GR_UNABLE";
}

void try_sqrt(gr_ctx_t ctx, const char * s)
{
    gr_ptr x, y;
    int status;

    GR_TMP_INIT2(x, y, ctx);
    GR_MUST_SUCCEED(gr_set_str(x, s, ctx));
    status = gr_sqrt(y, x, ctx);
    flint_printf("sqrt(%s) in %{gr_ctx}: %s", s, ctx, status_str(status));
    if (status == GR_SUCCESS)
        flint_printf(", %{gr}", y, ctx);
    flint_printf("\n");
    GR_TMP_CLEAR2(x, y, ctx);
}

int main(void)
{
    gr_ctx_t ZZ, RR, CC, QQbar, ZZ7;
    gr_ptr a, b;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_real_arb(RR, 53);
    gr_ctx_init_complex_acb(CC, 53);
    gr_ctx_init_real_qqbar(QQbar);
    gr_ctx_init_nmod(ZZ7, 7);

    try_sqrt(ZZ, "2");
    try_sqrt(QQbar, "2");
    try_sqrt(RR, "-1");
    try_sqrt(CC, "-1");
    try_sqrt(ZZ7, "2");
    try_sqrt(ZZ7, "3");
    try_sqrt(RR, "[0 +/- 0.1]");

    GR_TMP_INIT2(a, b, ZZ7);
    GR_MUST_SUCCEED(gr_set_ui(a, 2, ZZ7));
    flint_printf("exp(2) in Z/7Z: %s\n", status_str(gr_exp(b, a, ZZ7)));
    GR_TMP_CLEAR2(a, b, ZZ7);

    GR_TMP_INIT2(a, b, RR);
    GR_MUST_SUCCEED(gr_set_str(a, "1/3", RR));
    GR_MUST_SUCCEED(gr_mul_ui(b, a, 3, RR));
    flint_printf("3 * (1/3) == 1 in RR: ");
    truth_println(gr_is_one(b, RR));
    GR_MUST_SUCCEED(gr_set_str(b, "1/2", RR));
    flint_printf("1/3 == 1/2 in RR: ");
    truth_println(gr_equal(a, b, RR));
    GR_TMP_CLEAR2(a, b, RR);

    gr_ctx_clear(ZZ);
    gr_ctx_clear(RR);
    gr_ctx_clear(CC);
    gr_ctx_clear(QQbar);
    gr_ctx_clear(ZZ7);
    return 0;
}
```
```output
sqrt(2) in Integer ring (fmpz): GR_DOMAIN
sqrt(2) in Real algebraic numbers (qqbar): GR_SUCCESS, Root a = 1.41421 of a^2-2
sqrt(-1) in Real numbers (arb, prec = 53): GR_DOMAIN
sqrt(-1) in Complex numbers (acb, prec = 53): GR_SUCCESS, 1.000000000000000*I
sqrt(2) in Integers mod 7 (_gr_nmod): GR_SUCCESS, 3
sqrt(3) in Integers mod 7 (_gr_nmod): GR_DOMAIN
sqrt([0 +/- 0.1]) in Real numbers (arb, prec = 53): GR_UNABLE
exp(2) in Z/7Z: GR_UNABLE
3 * (1/3) == 1 in RR: T_UNKNOWN
1/3 == 1/2 in RR: T_FALSE
```

The square root of a ball that contains both negative and positive
numbers might or might not exist, hence `GR_UNABLE`. The exponential
function is not defined on $\Z/7\Z$ by any implementation, which is also
`GR_UNABLE`: the code says "I cannot", not "this has no value". Over
balls, $3 \cdot (1/3)$ is a ball that contains 1 but also other numbers,
so it is not known to equal 1, while $1/3$ and $1/2$ are separated and
known to be different. The algorithm page [](algorithms/generic/errors)
explains how algorithms use the distinction between the two failure
codes, and [](algorithms/generic/inexact) how generic algorithms work
with inexact elements.

## Vectors {#vectors}

Arrays of elements are stored contiguously; `GR_ENTRY(v, i, sz)` with
`sz = ctx->sizeof_elem` is the address of entry $i$. The functions of
[`gr_vec`](gr_vec) with a leading underscore act on such arrays with an
explicit length: elementwise arithmetic, dot products, sums and
products, conversions. The type `gr_vec_t` is a resizable vector that
manages its own memory, like `fmpz_poly_t`. Vector operations are not
just loops for convenience: rings override them with faster code (for
example dot products with delayed modular reduction), so a generic
algorithm should use them whenever it has a loop over elements.

```c example name=gr_vectors
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_vec.h>

int main(void)
{
    gr_ctx_t ZZ, ZZ7, QQ;
    gr_ptr a, b, s;
    gr_vec_t v;
    slong i, n = 5;
    int status = GR_SUCCESS;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_nmod(ZZ7, 7);
    gr_ctx_init_fmpq(QQ);

    /* a raw array of n integers: a[i] = i^2 */
    a = gr_heap_init_vec(n, ZZ);
    for (i = 0; i < n; i++)
        status |= gr_set_si(GR_ENTRY(a, i, ZZ->sizeof_elem), i * i, ZZ);

    /* convert to Z/7Z elementwise and take a dot product */
    b = gr_heap_init_vec(n, ZZ7);
    s = gr_heap_init(ZZ7);
    status |= _gr_vec_set_other(b, a, ZZ, n, ZZ7);
    status |= _gr_vec_dot(s, NULL, 0, b, b, n, ZZ7);
    flint_printf("a = %{gr*}\n", a, n, ZZ);
    flint_printf("b = %{gr*}, b.b = %{gr}\n", b, n, ZZ7, s, ZZ7);
    gr_heap_clear(s, ZZ7);

    /* a gr_vec_t manages its own length and memory */
    gr_vec_init(v, 0, QQ);
    status |= gr_vec_set_str(v, "[1/2, 1/3, -1/6]", 1, QQ);
    s = gr_heap_init(QQ);
    status |= gr_set_ui(s, 7, QQ);
    status |= gr_vec_append(v, s, QQ);
    status |= _gr_vec_sum(s, v->entries, gr_vec_length(v, QQ), QQ);
    flint_printf("v = %{gr*}, sum = %{gr}\n",
        v->entries, gr_vec_length(v, QQ), QQ, s, QQ);
    gr_heap_clear(s, QQ);
    gr_vec_clear(v, QQ);

    gr_heap_clear_vec(a, n, ZZ);
    gr_heap_clear_vec(b, n, ZZ7);
    gr_ctx_clear(ZZ);
    gr_ctx_clear(ZZ7);
    gr_ctx_clear(QQ);
    flint_printf("status %d\n", status);
    return 0;
}
```
```output
a = [0, 1, 4, 9, 16]
b = [0, 1, 4, 2, 2], b.b = 4
v = [1/2, 1/3, -1/6, 7], sum = 23/3
status 0
```

The second argument of `_gr_vec_dot` is an optional initial value
(`NULL` for none) and the third selects subtraction; see the function
page [](gr_vec/dot).

## Printing and parsing {#printing}

`gr_print` and `gr_println` print an element, `gr_get_str` returns a
string to be freed with `flint_free`, and `flint_printf` understands
`%{gr}` (element and context), `%{gr*}` (array, length and context),
`%{gr_ctx}` and `%{truth}` ([](flint/printf)). Real and complex balls
are printed with as many digits as are correct, `gr_get_str_n` rounds
to a given number of digits. In the other direction, `gr_set_str`
evaluates an expression with the ring's own arithmetic, so the same
string gives different values in different rings: `"1/3"` is a fraction
in $\mathbb{Q}$, $5$ in $\Z/7\Z$, a ball in $\mathbb{R}$, and `GR_DOMAIN` in
$\Z$. Polynomial and series rings recognise the names of their
generators (`x` by default, changed with `gr_ctx_set_gen_name`), and
ball notation `[m +/- r]` is accepted wherever balls are elements.

## Conversions between rings {#conversions}

`gr_set_other(res, x, x_ctx, ctx)` converts an element of one ring to
another. It succeeds when the conversion has a canonical meaning and is
implemented: $\Z$ to any ring, $\mathbb{Q}$ to rings in which the denominator is
invertible, $\Z[x]$ to $(\Z/7\Z)[x]$ as in the example above, a real
ball to a ball of another precision or to a complex ball. It returns
`GR_DOMAIN` when the element has no image (the rational $1/7$ in
$\Z/7\Z$, the complex number $i$ in $\mathbb{R}$) and `GR_UNABLE` when
no conversion is implemented: real balls, even exact ones, are not
converted to $\mathbb{Q}$ or $\Z$, and residues in $\Z/7\Z$ are not lifted
to $\Z$ (use `gr_get_fmpz` or the element type for that). Arithmetic between
elements of different rings is done by converting first, or with the
`_other` functions such as `gr_add_other` and `gr_mul_other`, which
convert one operand to the ring of the result.

## Where next {#next}

* [](gr/division): the different division functions and their failure
  modes.
* [](gr/implementing): how to define a new ring and test it.
* [`gr_poly`](gr_poly), [`gr_mat`](gr_mat), [`gr_series`](gr_series)
  and [`gr_mpoly`](gr_mpoly): generic polynomials, matrices, power series
  and multivariate polynomials.
* [](algorithms/generic/overview): the design of the framework, and
  [@Joh2025gr] for a detailed account.
