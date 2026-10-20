---
title: Zero tests and equality
functions: ca_check_is_zero ca_is_zero_check_fast ca_check_equal
---

Deciding whether a number is zero is the central difficulty of exact
real and complex arithmetic: it decides equality ($x = y$ if and only
if $x - y = 0$), division ($1/x$ is a number only if $x \ne 0$), the
degree of a polynomial, the rank of a matrix and the sign of a real
number. For the numbers that `ca_t` can represent the question is not
known to be decidable in general (see
[](algorithms/exact/zero-testing#decidability)), so
`ca_check_is_zero` returns a `truth_t`: `T_TRUE` and `T_FALSE` are
mathematically proved, and `T_UNKNOWN` means that the implemented
methods failed to decide.

::: summary algorithms/exact/zero-testing :::

## How the answer is obtained {#strategy}

Proving that a number is nonzero is easy if it is: evaluating it in
ball arithmetic at sufficiently high precision gives a ball that
excludes zero. Proving that a number is zero requires an exact argument.
`ca_check_is_zero` combines these as follows.

1. The representation is inspected (`ca_is_zero_check_fast`). Special
   values are not zero, and an element of $\mathbb{Q}$ or of a field
   $\mathbb{Q}(a)$ with a single algebraic generator is zero if and only
   if its canonical representation is zero.
2. In a field with several generators, the denominator is dropped and
   the numerator is evaluated in ball arithmetic at 64, 128, ... bits up
   to `CA_OPT_PREC_LIMIT`. A ball that excludes zero proves $x \ne 0$.
3. At the first precision, if all generators are algebraic (square
   roots, roots and complex parts of algebraic numbers), $x$ is
   evaluated exactly as an algebraic number with `ca_get_qqbar`, which
   decides the question if the degrees stay within
   `CA_OPT_QQBAR_DEG_LIMIT`.
4. If the numerical evaluation is still inconclusive at the maximum
   precision, the generators are rewritten recursively in terms of
   exponentials and logarithms (trigonometric functions become complex
   exponentials, complex parts become conjugates, algebraic numbers
   become elements of cyclotomic fields where possible), which often
   exposes relations that the original field does not record; a
   representation that becomes zero proves $x = 0$.
5. Finally the numerator is factored as a polynomial in the
   generators; if it has at least two nonconstant factors, each factor
   is tested by steps 1 to 4, and if one of them is proved zero, so is
   $x$.

Relations between the generators that are known when a field is
created, such as $(\sqrt{2})^2 = 2$, $e^{a} e^{-a} = 1$ or $\log 8 = 3
\log 2$ (see [](algorithms/exact/field-towers)), are applied by every
arithmetic operation, so many zero results are recognised in step 1
before any test is made.

`ca_check_equal` compares representations first, decides comparisons
with rational numbers and within one number field directly, tries to
separate $x$ and $y$ numerically at 64 bits, and otherwise tests $x - y$
with `ca_check_is_zero`. The predicates `ca_check_is_one`,
`ca_check_is_integer`, `ca_check_is_real` and the comparisons
`ca_check_lt` and so on are built in the same way.

## Example {#example}

The program below tests four numbers. The first is a denesting identity
of algebraic numbers, proved by exact algebraic evaluation. The second
is Ramanujan's near-integer $e^{\pi \sqrt{163}} \approx 640320^3 + 744 -
7.5 \cdot 10^{-13}$, which numerical evaluation proves nonzero. The third
is zero because $(\sqrt{2} + 1)(\sqrt{2} - 1) = 1$; the integer relation
between the two logarithms is found and proved when their common field
is created, so the sum simplifies to zero. The last is zero by the
reflection formula $\Gamma(s) \Gamma(1 - s) = \pi / \sin(\pi s)$, which
the implementation does not know, so the answer is `T_UNKNOWN`.

```c example name=ca_check_is_zero
#include <flint/calcium.h>
#include <flint/ca.h>

static void report(const char * descr, const ca_t x, ca_ctx_t ctx)
{
    flint_printf("%-42s ", descr);
    truth_print(ca_check_is_zero(x, ctx));
    flint_printf("\n");
}

int main(void)
{
    ca_ctx_t ctx;
    ca_t x, y, z;

    ca_ctx_init(ctx);
    ca_init(x, ctx);
    ca_init(y, ctx);
    ca_init(z, ctx);

    /* sqrt(5 + 2 sqrt(6)) - sqrt(2) - sqrt(3) */
    ca_sqrt_ui(x, 6, ctx);
    ca_mul_ui(x, x, 2, ctx);
    ca_add_ui(x, x, 5, ctx);
    ca_sqrt(x, x, ctx);
    ca_sqrt_ui(y, 2, ctx);
    ca_sub(x, x, y, ctx);
    ca_sqrt_ui(y, 3, ctx);
    ca_sub(x, x, y, ctx);
    report("sqrt(5+2*sqrt(6)) - sqrt(2) - sqrt(3)", x, ctx);

    /* exp(pi sqrt(163)) - (640320^3 + 744) */
    ca_pi(x, ctx);
    ca_sqrt_ui(y, 163, ctx);
    ca_mul(x, x, y, ctx);
    ca_exp(x, x, ctx);
    ca_set_ui(y, 640320, ctx);
    ca_pow_ui(y, y, 3, ctx);
    ca_add_ui(y, y, 744, ctx);
    ca_sub(x, x, y, ctx);
    report("exp(pi*sqrt(163)) - (640320^3 + 744)", x, ctx);

    /* log(sqrt(2) + 1) + log(sqrt(2) - 1) */
    ca_sqrt_ui(y, 2, ctx);
    ca_add_ui(x, y, 1, ctx);
    ca_log(x, x, ctx);
    ca_sub_ui(y, y, 1, ctx);
    ca_log(y, y, ctx);
    ca_add(x, x, y, ctx);
    report("log(sqrt(2) + 1) + log(sqrt(2) - 1)", x, ctx);

    /* Gamma(1/3) Gamma(2/3) - 2 pi / sqrt(3) */
    ca_set_ui(x, 1, ctx);
    ca_div_ui(x, x, 3, ctx);
    ca_gamma(y, x, ctx);
    ca_ui_sub(x, 1, x, ctx);
    ca_gamma(z, x, ctx);
    ca_mul(y, y, z, ctx);
    ca_pi(x, ctx);
    ca_mul_ui(x, x, 2, ctx);
    ca_sqrt_ui(z, 3, ctx);
    ca_div(x, x, z, ctx);
    ca_sub(x, y, x, ctx);
    report("Gamma(1/3) Gamma(2/3) - 2 pi / sqrt(3)", x, ctx);

    ca_clear(x, ctx);
    ca_clear(y, ctx);
    ca_clear(z, ctx);
    ca_ctx_clear(ctx);
    return 0;
}
```
```output
sqrt(5+2*sqrt(6)) - sqrt(2) - sqrt(3)      T_TRUE
exp(pi*sqrt(163)) - (640320^3 + 744)       T_FALSE
log(sqrt(2) + 1) + log(sqrt(2) - 1)        T_TRUE
Gamma(1/3) Gamma(2/3) - 2 pi / sqrt(3)     T_UNKNOWN
```

## Dealing with T_UNKNOWN {#unknown}

`T_UNKNOWN` is returned when $x$ cannot be separated from zero by
numerical evaluation up to `CA_OPT_PREC_LIMIT` bits and no exact
argument applies. In practice this almost always means that $x = 0$ (or
that $x$ is nonzero but extremely small) and that the identity involves
transcendental numbers related in a way the implementation does not
know, such as the reflection formula above. `T_UNKNOWN` is not an
error but a legitimate answer, and code that uses exact numbers must
handle it.

The following program proves Machin's formula, and then tests the
nonzero number $\pi - p/2^{5000}$, where $p/2^{5000}$ is $\pi$ truncated
to 5000 bits. The difference is about $2^{-5000}$, too small to be
separated from zero at the default limit of 4096 bits; raising the limit
with `ca_ctx_set_option` proves it nonzero.

```c example name=ca_check_is_zero_limits
#include <flint/arb.h>
#include <flint/calcium.h>
#include <flint/ca.h>

int main(void)
{
    ca_ctx_t ctx;
    ca_t x, y;
    arb_t t;
    fmpz_t p, q;
    truth_t r;

    ca_ctx_init(ctx);
    ca_init(x, ctx);
    ca_init(y, ctx);
    arb_init(t);
    fmpz_init(p); fmpz_init(q);
    fmpz_one_2exp(q, 5000);

    /* Machin's formula: 16 atan(1/5) - 4 atan(1/239) - pi */
    ca_set_ui(x, 5, ctx); ca_inv(x, x, ctx); ca_atan(x, x, ctx);
    ca_set_ui(y, 239, ctx); ca_inv(y, y, ctx); ca_atan(y, y, ctx);
    ca_mul_ui(x, x, 16, ctx);
    ca_mul_ui(y, y, 4, ctx);
    ca_sub(x, x, y, ctx);
    ca_pi(y, ctx);
    ca_sub(x, x, y, ctx);
    flint_printf("Machin: "); truth_print(ca_check_is_zero(x, ctx)); flint_printf("\n");

    /* pi - p / 2^5000 with p = floor(pi 2^5000): about 2^-5000 */
    arb_const_pi(t, 5100);
    arb_mul_2exp_si(t, t, 5000);
    arf_get_fmpz(p, arb_midref(t), ARF_RND_FLOOR);
    ca_set_fmpz(y, p, ctx);
    ca_div_fmpz(y, y, q, ctx);
    ca_pi(x, ctx);
    ca_sub(x, x, y, ctx);

    r = ca_check_is_zero(x, ctx);
    flint_printf("pi - p/2^5000, limit 4096: "); truth_print(r); flint_printf("\n");

    ca_ctx_set_option(ctx, CA_OPT_PREC_LIMIT, 8192);
    r = ca_check_is_zero(x, ctx);
    flint_printf("pi - p/2^5000, limit 8192: "); truth_print(r); flint_printf("\n");

    fmpz_clear(p); fmpz_clear(q);
    arb_clear(t);
    ca_clear(x, ctx);
    ca_clear(y, ctx);
    ca_ctx_clear(ctx);
    return 0;
}
```

```output
Machin: T_TRUE
pi - p/2^5000, limit 4096: T_UNKNOWN
pi - p/2^5000, limit 8192: T_FALSE
```

Raising `CA_OPT_PREC_LIMIT` turns `T_UNKNOWN` into `T_FALSE` for tiny
nonzero numbers; it can lead to `T_TRUE` only indirectly, since the
exact algebraic evaluation limits bit sizes to ten times
`CA_OPT_PREC_LIMIT`. In general, ways to proceed are:

- if the number is expected to be nonzero but very small, raise
  `CA_OPT_PREC_LIMIT`;
- if the number is algebraic of high degree, raise
  `CA_OPT_QQBAR_DEG_LIMIT`, or compute with [](qqbar) directly;
- if the number involves transcendental functions and is expected to be
  zero, rewrite the input so that the relation is visible: for example,
  express trigonometric functions through exponentials (the default
  `CA_TRIG_EXPONENTIAL` form does this), or express $\Gamma(2/3)$
  through $\Gamma(1/3)$ with the reflection formula before computing;
- treat the result as a numerical certificate if the application
  permits ("zero to 4096 bits"), keeping in mind that this is not a
  proof;
- report the result as undecided, which is what functions of [](ca),
  [](ca_poly) and [](ca_mat) do: an operation that depends on an
  undecided zero test returns *Unknown*, `T_UNKNOWN` or a failure flag.

The options are described in [](ca/options); the limits apply to the
whole context, so they can be set once before a computation. Since a
context caches fields and enclosures, repeating a test in the same
context is cheap.
