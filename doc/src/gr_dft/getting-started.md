---
title: Getting started
summary: The conventions of the generic rings module that gr_dft relies on, which rings can be used, and a first transform over complex balls and over Z/pZ.
order: 1
---

`gr_dft` is written on top of `gr`, FLINT's interface to generic rings, so
one implementation serves complex balls, exact algebraic numbers,
$\Z/m\Z$, finite fields and fixed-point numbers.

## Generic rings in brief {#gr}

**Contexts.** A ring is described by a context of type `gr_ctx_t`,
initialised by a constructor such as `gr_ctx_init_nmod(ctx, m)` for
$\Z/m\Z$ with word-size $m$, `gr_ctx_init_fmpz_mod(ctx, m)` for an
`fmpz` modulus, `gr_ctx_init_fq(ctx, p, d, "a")` for a finite field,
`gr_ctx_init_complex_acb(ctx, prec)` and `gr_ctx_init_real_arb(ctx, prec)`
for complex and real balls with working precision `prec` bits,
`gr_ctx_init_complex_qqbar(ctx)` for exact algebraic numbers, or
`gr_dft_ctx_init_nfixed_complex(ctx, nlimbs)` from this module. Every
context is freed with `gr_ctx_clear`. Functions of `gr` take the context
as their last argument.

**Elements.** An element is an untyped block of `ctx->sizeof_elem` bytes,
addressed by a `gr_ptr` (the output) or a `gr_srcptr` (an input, `const`).
A vector is a contiguous array of elements, and `GR_ENTRY(vec, i, sz)`
with `sz = ctx->sizeof_elem` is the address of entry $i$. Vectors are
allocated and initialised with `gr_heap_init_vec(len, ctx)` and freed
with `gr_heap_clear_vec(vec, len, ctx)`; single elements with
`gr_heap_init` and `gr_heap_clear`, or on the stack with `GR_TMP_INIT`
and `GR_TMP_CLEAR`. Elements of a specific type can be used directly: an
element of an `acb` context is an `acb_struct`, so an `acb_ptr` can be
passed where a `gr_ptr` is expected, and an element of an `nmod` context
is a `ulong`.

**Status codes.** Almost every function returns an `int` status:

| Value | Meaning | Example |
|---|---|---|
| `GR_SUCCESS` (0) | the result is correct | |
| `GR_DOMAIN` (1) | the result is not defined in this ring | division by $n$ modulo $m$ with $\gcd(n, m) > 1$; $e^{-2\pi i/3}$ in $\mathbb{R}$ |
| `GR_UNABLE` (2) | this implementation cannot compute the result | no method to compute $e^{-2\pi i/n}$ over $\Z/m\Z$; overflow in fixed point |

Statuses are combined with `|`, so `GR_DOMAIN | GR_UNABLE` (3) can occur,
and code typically accumulates them with `status |= f(...)`. The functions
of this module that can fail are declared `WARN_UNUSED_RESULT`, so the
compiler warns when a status is ignored. When failure is impossible or
fatal, wrap a call in `GR_MUST_SUCCEED`, which aborts with an error
message; `GR_IGNORE` discards a status explicitly. After a failure the
output values are unspecified.

**Predicates.** Comparisons such as `gr_equal` and `gr_is_one` return a
`truth_t`: `T_TRUE`, `T_FALSE` or `T_UNKNOWN`, the last when an inexact
ring cannot decide. `T_TRUE` is 0, so the result must be
compared with `T_TRUE` rather than used as a C boolean.

## Which rings can be used {#rings}

| Ring | Context | Roots of unity | Remarks |
|---|---|---|---|
| complex balls | `gr_ctx_init_complex_acb` | default | use `gr_dft_precomp_init_karatsuba` with `gr_ctx_init_real_arb`, or see [](gr_dft/acb) |
| algebraic numbers | `gr_ctx_init_complex_qqbar`, `gr_ctx_init_complex_ca` | default | exact but slow |
| fixed point | `gr_dft_ctx_init_nfixed_complex` | default | inputs must be scaled; see [](gr_dft/nfixed) |
| $\Z/m\Z$, $\mathbb{F}_q$ | `gr_ctx_init_nmod`, `gr_ctx_init_fmpz_mod`, `gr_ctx_init_fq` | supplied by the caller | $n \mid q - 1$ for a field $\mathbb{F}_q$ |
| $\Z$, $\mathbb{Q}$, $\mathbb{R}$ | `gr_ctx_init_fmpz`, `gr_ctx_init_fmpq`, `gr_ctx_init_real_arb` | only $n \le 2$ | |

Over rings without default roots, the one-shot functions `gr_dft`,
`gr_dft_inverse` and the product and Dirichlet transforms return
`GR_UNABLE` or `GR_DOMAIN`, and a plan must be made with
`gr_dft_precomp_init_root` (or `gr_dft_prod_precomp_init_root`). [](gr_dft/plans) explains how to choose the root and what the
constructors check.

## A first transform {#example}

The program computes the transform of $(1, 2, 3, 4)$ over complex balls
with the default root $-i$, and then over $\Z/17\Z$, where `gr_dft` fails
because there is no default root and the root 4, of order 4, is supplied
instead.

```c example name=gr_dft_first
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_vec.h>
#include <flint/gr_dft.h>

int main(void)
{
    gr_ctx_t CC, ZP;
    gr_ptr x, y, w;
    slong j, n = 4;
    int status;

    /* complex balls with 53-bit precision: the default root works */
    gr_ctx_init_complex_acb(CC, 53);
    x = gr_heap_init_vec(n, CC);
    y = gr_heap_init_vec(n, CC);
    for (j = 0; j < n; j++)
        GR_MUST_SUCCEED(gr_set_si(GR_ENTRY(x, j, CC->sizeof_elem), j + 1, CC));

    status = gr_dft(y, x, n, CC);
    flint_printf("status %d: ", status);
    _gr_vec_print(y, n, CC);
    flint_printf("\n");
    gr_heap_clear_vec(x, n, CC);
    gr_heap_clear_vec(y, n, CC);
    gr_ctx_clear(CC);

    /* Z/17Z: there is no default root, so supply w = 4 of order 4 */
    gr_ctx_init_nmod(ZP, 17);
    x = gr_heap_init_vec(n, ZP);
    y = gr_heap_init_vec(n, ZP);
    w = gr_heap_init(ZP);
    for (j = 0; j < n; j++)
        GR_MUST_SUCCEED(gr_set_si(GR_ENTRY(x, j, ZP->sizeof_elem), j + 1, ZP));

    status = gr_dft(y, x, n, ZP);
    flint_printf("gr_dft over Z/17Z: %s\n",
        status == GR_UNABLE ? "GR_UNABLE" : "?");

    {
        gr_dft_pre_t P;

        GR_MUST_SUCCEED(gr_set_ui(w, 4, ZP));
        status = gr_dft_precomp_init_root(P, w, n, GR_DFT_ALG_AUTO, 0, ZP);
        if (status == GR_SUCCESS)
        {
            status |= gr_dft_precomp(y, x, P, ZP);
            gr_dft_precomp_clear(P);
        }
        flint_printf("status %d: ", status);
        _gr_vec_print(y, n, ZP);
        flint_printf("\n");
    }

    gr_heap_clear(w, ZP);
    gr_heap_clear_vec(x, n, ZP);
    gr_heap_clear_vec(y, n, ZP);
    gr_ctx_clear(ZP);
    return 0;
}
```
```output
status 0: [10.00000000000000, (-2.000000000000000 + 2.000000000000000*I), -2.000000000000000, (-2.000000000000000 - 2.000000000000000*I)]
gr_dft over Z/17Z: GR_UNABLE
status 0: [10, 7, 15, 6]
```

Over $\Z/17\Z$, for example, $X_1 = 1 + 2 \cdot 4 + 3 \cdot 4^2 + 4 \cdot 4^3
= 313 \equiv 7$. The balls over $\mathbb{C}$ are exact here, because the
roots of unity of order 4 are computed exactly.

## Where next {#next}

- [](gr_dft/plans) and [](gr_dft/transforms): plans, roots, status codes
  and convolutions.
- [](gr_dft/acb): rigorous transforms of complex ball vectors.
- [](gr_dft/prod): multidimensional transforms and Dirichlet characters.
- [](gr_dft/performance): algorithms, complex mode and threads.
