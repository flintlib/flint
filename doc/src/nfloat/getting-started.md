---
title: Getting started with nfloat
summary: Creating nfloat contexts, declaring elements on the stack, status codes, directed rounding, and running generic algorithms such as matrix inversion over nfloat.
order: 1
---

An `nfloat` context fixes the precision, a whole number of words, and
the treatment of special values. Elements are plain arrays of words:
they can be declared on the stack with the `nfloatX_t` types (or
allocated with `GR_TMP_INIT` or `gr_heap_init` for a precision chosen at
run time) and need no clearing. Every function returns a status code,
`GR_SUCCESS` (0) when the result is valid; operations whose result is
not a finite real number, such as $\log(-1)$, division by zero, or
overflow, return `GR_DOMAIN` or `GR_UNABLE` unless the context was
created with flags that allow infinities and NaN.

## Basic use {#basic}

The program below computes $\sqrt 2$ and $e$ at 128 bits, shows the
status of an invalid operation, and brackets $1/3$ between the results
of two contexts with directed rounding. Since `nfloat_get_arf` is only
declared when `arf.h` is included, that header comes first.

```c example name=nfloat_basic
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/arf.h>
#include <flint/nfloat.h>

int main(void)
{
    gr_ctx_t ctx, lo, hi;
    nfloat128_t x, y;
    nfloat64_t a, b;
    arf_t s, t;
    int status = GR_SUCCESS;

    /* 128-bit floating-point numbers on the stack */
    status |= nfloat_ctx_init(ctx, 128, 0);
    gr_ctx_println(ctx);
    status |= nfloat_set_ui(x, 2, ctx);
    status |= nfloat_sqrt(x, x, ctx);
    status |= gr_println(x, ctx);
    status |= nfloat_one(y, ctx);
    status |= nfloat_exp(y, y, ctx);
    status |= gr_println(y, ctx);

    /* errors are reported through status codes */
    status |= nfloat_set_si(y, -1, ctx);
    flint_printf("log(-1): %s\n",
        nfloat_log(y, y, ctx) == GR_SUCCESS ? "success" : "failure");

    /* directed rounding brackets the exact quotient 1/3 */
    status |= nfloat_ctx_init(lo, 64, NFLOAT_RND_FLOOR);
    status |= nfloat_ctx_init(hi, 64, NFLOAT_RND_CEIL);
    status |= nfloat_set_ui(a, 1, lo);
    status |= nfloat_div_ui(a, a, 3, lo);
    status |= nfloat_set_ui(b, 1, hi);
    status |= nfloat_div_ui(b, b, 3, hi);
    arf_init(s);
    arf_init(t);
    status |= nfloat_get_arf(s, a, lo);
    status |= nfloat_get_arf(t, b, hi);
    arf_sub(t, t, s, ARF_PREC_EXACT, ARF_RND_DOWN);
    flint_printf("upper - lower = ");
    arf_printd(t, 5);
    flint_printf("\nstatus = %d\n", status);

    arf_clear(s);
    arf_clear(t);
    gr_ctx_clear(ctx);
    gr_ctx_clear(lo);
    gr_ctx_clear(hi);
    return 0;
}
```
```output
Floating-point numbers with prec = 128 (nfloat)
1.41421356237309504880168872420969807857
2.71828182845904523536028747135266249775
log(-1): failure
upper - lower = 2.7105e-20
status = 0
```

The two bounds for $1/3 \in [2^{-2}, 2^{-1})$ differ by
$2^{-65} \approx 2.7 \cdot 10^{-20}$, one unit in the last place at 64
bits.

## Generic algorithms {#generic}

Since an `nfloat` context is a generic ring, the polynomial, matrix and
power series algorithms of `gr_poly`, `gr_mat` and `gr_series` work over
it. The context overrides the operations that matter most for speed:
vector arithmetic and dot products (`_nfloat_vec_dot`), matrix
multiplication (`nfloat_mat_mul`), and the triangular solving and LU
factorisation on which `gr_mat_inv`, `gr_mat_det` and
`gr_mat_nonsingular_solve` are built.

The following program inverts the $10 \times 10$ Hilbert matrix at 192
bits. The exact inverse has integer entries; the last one is
$44914183600$. The result is approximate, with a relative error of
about $6 \cdot 10^{-47}$: the condition number of the matrix, about
$1.6 \cdot 10^{13}$, costs some 12 of the 57 digits of the working
precision. To obtain a
certified result, use `arb_mat_t`, or the `dfloat` ball contexts at up
to 212 bits.

```c example name=nfloat_hilbert_inverse
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_mat.h>
#include <flint/nfloat.h>

int main(void)
{
    gr_ctx_t ctx;
    gr_mat_t H, Hinv;
    slong i, j, n = 10, sz;
    int status = GR_SUCCESS;

    status |= nfloat_ctx_init(ctx, 192, 0);
    sz = ctx->sizeof_elem;

    gr_mat_init(H, n, n, ctx);
    gr_mat_init(Hinv, n, n, ctx);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            status |= gr_set_ui(GR_MAT_ENTRY(H, i, j, sz), i + j + 1, ctx);
            status |= gr_inv(GR_MAT_ENTRY(H, i, j, sz), GR_MAT_ENTRY(H, i, j, sz), ctx);
        }

    status |= gr_mat_inv(Hinv, H, ctx);
    flint_printf("Hinv[9][9] = ");
    status |= gr_println(GR_MAT_ENTRY(Hinv, 9, 9, sz), ctx);
    flint_printf("status = %d\n", status);

    gr_mat_clear(H, ctx);
    gr_mat_clear(Hinv, ctx);
    gr_ctx_clear(ctx);
    return 0;
}
```
```output
Hinv[9][9] = 44914183599.99999999999999999999999999999999999729212779167
status = 0
```

## Choosing the precision and the function precision {#precision}

The precision is rounded up to a multiple of 64 bits (on 64-bit
machines), and operations at one and two words have the fastest
kernels, with most of the arithmetic in registers. When only some of the
bits are needed for transcendental functions, for example when 53-bit
results suffice but the arithmetic needs 64 bits, a smaller function
precision set with `nfloat_ctx_set_func_prec` makes the elementary
functions cheaper without changing the arithmetic.

The complex numbers `nfloat_complex` are created with
`nfloat_complex_ctx_init` and work in the same way, with the real and
imaginary parts stored consecutively.
