---
title: Context options
functions: ca_ctx_get_option CA_OPT_VERBOSE CA_OPT_MPOLY_ORD CA_OPT_PREC_LIMIT CA_OPT_QQBAR_DEG_LIMIT CA_OPT_USE_GROEBNER CA_OPT_VIETA_LIMIT CA_OPT_TRIG_FORM CA_OPT_NUM_OPTIONS
---

Every `ca_ctx_t` carries an array of options that control how hard the
functions of [](ca) try to simplify and to decide predicates, how
trigonometric functions are represented, and how numbers are printed.
They are read with `ca_ctx_get_option` and set with `ca_ctx_set_option`
(or directly as `ctx->options[i]`), using the indices `CA_OPT_...`.
Since fields, relations and numerical enclosures are cached in the
context, options that affect evaluation should be set right after
`ca_ctx_init`, before any number is created; a later change does not
apply to objects that are already cached, and two numbers computed
before and after a change may be represented differently. Printing
options can be changed at any time.

## Evaluation limits {#limits}

The limits trade the ability to simplify and to decide predicates
against running time. None of them affects correctness: a result that
is returned is always correct, but with smaller limits more predicates
return `T_UNKNOWN` and more operations return *Unknown* or fail.

| Option | Default | Effect |
|---|---|---|
| `CA_OPT_PREC_LIMIT` | 4096 | Maximum precision in bits for numerical evaluation, for example in zero tests and comparisons; also bounds the size of exact powers and of qqbar computations (about ten times this many bits) |
| `CA_OPT_LOW_PREC` | 64 | Precision for quick numerical checks, such as the numerical coprimality test in `ca_poly_gcd` |
| `CA_OPT_LLL_PREC` | 128 | Precision for finding integer relations between logarithms and between exponentials and powers with LLL when a field is created |
| `CA_OPT_QQBAR_DEG_LIMIT` | 120 | Maximum degree of algebraic numbers computed internally, for example by `ca_get_qqbar`, in simplifications of square roots and complex parts, and for roots of unity $e^{\pi i p/q}$ ($q$ at most this value) |
| `CA_OPT_POW_LIMIT` | 20 | Largest exponent $n$ for which $x^n$ is expanded by multiplication in generic fields; larger powers become extension numbers $x^n$ |
| `CA_OPT_SMOOTH_LIMIT` | 32 | Size in bits of the prime factors found by smooth factorisation, used when simplifying square roots |
| `CA_OPT_USE_GROEBNER` | 1 | Whether to replace the ideal of relations of each new field by its Gröbner basis |
| `CA_OPT_GROEBNER_LENGTH_LIMIT` | 100 | Maximum number of polynomials in Buchberger's algorithm |
| `CA_OPT_GROEBNER_POLY_LENGTH_LIMIT` | 1000 | Maximum number of terms of a polynomial in Buchberger's algorithm |
| `CA_OPT_GROEBNER_POLY_BITS_LIMIT` | 10000 | Maximum coefficient size in bits in Buchberger's algorithm |
| `CA_OPT_VIETA_LIMIT` | 6 | Maximum degree $n$ of an algebraic number for which Vieta's formulas are added when all its conjugates are generators of one field; the formulas have $O(2^n)$ terms |

When Buchberger's algorithm exceeds a limit, the field keeps the
relations found so far, which are still valid but may reduce elements
less completely. The algorithms that use these options are described in
[](algorithms/exact/field-towers) and [](algorithms/exact/zero-testing).

## Representation {#representation}

`CA_OPT_MPOLY_ORD` selects the monomial order (`ORD_LEX`, the default,
`ORD_DEGLEX` or `ORD_DEGREVLEX`) of the multivariate polynomials that
represent elements of generic fields. The lexicographic order, combined
with the order of the generators, makes reduction modulo the ideal of
relations eliminate the most complicated generators. The option must be
set before any computation, since the polynomial contexts are created
on first use.

`CA_OPT_TRIG_FORM` selects how `ca_sin_cos`, `ca_tan` and the inverse
trigonometric functions represent their results:

- `CA_TRIG_EXPONENTIAL` (default): in terms of $e^{ix}$ and complex
  logarithms. All trigonometric values of one argument then live in one
  field $\mathbb{Q}(i, e^{ix})$, identities such as $\sin^2 x + \cos^2 x
  = 1$ hold by construction, and values at rational multiples of $\pi$
  become elements of cyclotomic fields. The price is that real numbers
  are written with complex generators, so that their numerical
  enclosures may have a tiny nonzero imaginary part.
- `CA_TRIG_DIRECT`: as the functions themselves, $\sin x$, $\cos x$,
  $\tan x$, $\operatorname{atan} x$, after symmetry reductions. Real
  inputs give real generators, but fewer relations are known.
- `CA_TRIG_SINE_COSINE`: like the direct form, with $\tan x$ written as
  $\sin x / \cos x$.
- `CA_TRIG_TANGENT`: sines and cosines in terms of $t = \tan(x/2)$, as
  $2t/(1 + t^2)$ and $(1 - t^2)/(1 + t^2)$.

## Printing and debugging {#printing}

`CA_OPT_PRINT_FLAGS` is the style used by `ca_print`, `ca_get_str` and
the printing functions of the vector, polynomial and matrix modules; it
combines the flags `CA_PRINT_N`, `CA_PRINT_REPR`, `CA_PRINT_FIELD` and
a number of digits $d$ given as $d \cdot$ `CA_PRINT_DIGITS` (see
[](ca#printing)). `CA_OPT_VERBOSE`, if nonzero, makes some internal
algorithms print what they do, which helps to understand why a
predicate returns `T_UNKNOWN`.

## Example {#example}

The program below prints $\sin(1)$ in the default exponential form,
then with 30 digits and without the symbolic part, and finally in a
second context that uses the direct form. The tiny imaginary part shown
in the first two lines comes from evaluating a real number through
complex exponentials; it is an artefact of the enclosure, not of the
exact value.

```c example name=ca_ctx_options
#include <flint/ca.h>

int main(void)
{
    ca_ctx_t ctx, ctx2;
    ca_t x, s, t;

    ca_ctx_init(ctx);
    ca_init(x, ctx);
    ca_init(s, ctx);

    /* the default trigonometric form uses complex exponentials */
    ca_set_ui(x, 1, ctx);
    ca_sin(s, x, ctx);
    flint_printf("sin(1) = "); ca_print(s, ctx); flint_printf("\n");

    /* print 30 digits, without the symbolic representation */
    ca_ctx_set_option(ctx, CA_OPT_PRINT_FLAGS, CA_PRINT_N | (CA_PRINT_DIGITS * 30));
    flint_printf("sin(1) = "); ca_print(s, ctx); flint_printf("\n");

    /* a second context with direct trigonometric functions */
    ca_ctx_init(ctx2);
    ca_ctx_set_option(ctx2, CA_OPT_TRIG_FORM, CA_TRIG_DIRECT);
    ca_init(t, ctx2);
    ca_transfer(t, ctx2, x, ctx);
    ca_sin(t, t, ctx2);
    flint_printf("sin(1) = "); ca_print(t, ctx2); flint_printf("\n");

    ca_clear(t, ctx2);
    ca_ctx_clear(ctx2);
    ca_clear(x, ctx);
    ca_clear(s, ctx);
    ca_ctx_clear(ctx);
    return 0;
}
```
```output
sin(1) = 0.841471 + 0e-24*I {(-a^2*b+b)/(2*a) where a = 0.540302 + 0.841471*I [Exp(1.00000*I {b})], b = I [b^2+1=0]}
sin(1) = 0.841470984807896506652502321630 - 0e-48*I
sin(1) = 0.841471 {a where a = 0.841471 [Sin(1)]}
```

Each number must be cleared with the context it belongs to, before that
context is cleared; `ca_transfer` moves values between contexts.
