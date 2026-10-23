---
title: Using dfloat
summary: The direct interface and the generic rings, with examples of plain expansions, exact ball arithmetic, status codes and a linear system solved with triple-double balls.
order: 1
---

The `dfloat` formats can be used in two ways. The functions with a
format prefix (`d2_add`, `d3b_exp`, ...) operate on stack-allocated
structs and need no initialisation or clearing: a `d2_t` is simply two
doubles. Alternatively, `gr_ctx_init_dfloat` creates a generic ring
whose elements have the same layout, so that every algorithm of the
generic modules (polynomials, matrices, power series) runs on the
formats with the fast vector kernels of this module.

## The direct interface {#direct}

A plain expansion is an unevaluated sum of doubles. Values that need
more than 53 bits are represented exactly as long as they fit in $53N$
bits, and the arithmetic is accurate to a few ulps at that precision.
For the balls, every function either returns nothing (it cannot fail)
or a status code: `GR_SUCCESS` when the result is a valid enclosure,
`GR_DOMAIN` when the input is certainly outside the domain of the
function, and `GR_UNABLE` when the input balls straddle the boundary of
the domain (see [](dfloat/balls)).

```c example name=dfloat_basic
#include <flint/flint.h>
#include <flint/dfloat.h>

int main(void)
{
    d2_t x, y;
    d2b_t a, b;
    int status;

    if (!dfloat_is_supported())
        return 0;

    /* 1 + 2^-60 is exact in double-double */
    d2_set_dd(x, 1.0, 0x1p-60);
    d2_set_d(y, 1.0);
    d2_sub(x, x, y);
    flint_printf("(1 + 2^-60) - 1 = %.17g\n", d2_get_d(x));

    /* exp(1) to about 106 bits */
    d2_one(x);
    d2_exp(x, x);
    d2_print(x); flint_printf("\n");

    /* balls: exact operations on exact balls stay exact */
    d2b_set_d(a, 3.0);
    d2b_mul(b, a, a);
    flint_printf("3*3 exact: %d\n", d2b_is_exact(b));

    /* 1/3 is rounded; 3 * (1/3) is a ball containing 1 */
    status = d2b_inv(b, a);
    d2b_mul(b, b, a);
    d2b_print(b);
    flint_printf("  status %d, contains 1: %d\n", status, d2b_contains_d(b, 1.0));

    /* division by a ball containing zero is undecided */
    d2b_set_d_rad(a, 0.0, 1.0);
    d2b_one(b);
    status = d2b_div(b, b, a);
    flint_printf("1 / [0 +/- 1]: status %d (GR_UNABLE = %d)\n", status, GR_UNABLE);
    return 0;
}
```
```output
(1 + 2^-60) - 1 = 8.6736173798840355e-19
2.7182818284590452353602874713527
3*3 exact: 1
[1.0000000000000000000000000000000 +/- 1.24e-32]  status 0, contains 1: 1
1 / [0 +/- 1]: status 2 (GR_UNABLE = 2)
```

## Generic rings {#rings}

`gr_ctx_init_dfloat(ctx, n, flags)` selects the number of components
and, with the flags, balls (`DFLOAT_BALL`), complex numbers
(`DFLOAT_COMPLEX`), canonical results (`DFLOAT_STRONG`) or cheaper
error bounds for ball products (`DFLOAT_FAST`). It returns `GR_UNABLE`
on targets where the module does not work (see `dfloat_is_supported`),
so that portable code can fall back to `arb_t`. The rings implement the
arithmetic and the elementary functions; special functions such as the
gamma function return `GR_UNABLE`.

The following program evaluates $e^\pi$ in `d3b` balls and solves a
$6 \times 6$ Hilbert system $Hx = (1, \ldots, 1)^T$, whose exact
solution is an integer vector. The radii, about $10^{-34}$, show that
the triple-double balls carry about 113 bits through Gaussian
elimination on a matrix with condition number $1.5 \cdot 10^7$.

```c example name=dfloat_gr_hilbert
#include <flint/flint.h>
#include <flint/arb.h>
#include <flint/gr.h>
#include <flint/gr_mat.h>
#include <flint/gr_special.h>
#include <flint/dfloat.h>

int main(void)
{
    gr_ctx_t ctx, RR;
    gr_mat_t H, X, B;
    gr_ptr x;
    arb_t y;
    slong i, j, n = 6, sz;
    int status = GR_SUCCESS;

    if (gr_ctx_init_dfloat(ctx, 3, DFLOAT_BALL) != GR_SUCCESS)
        return 0;   /* dfloat is not supported on this target */
    sz = ctx->sizeof_elem;

    gr_ctx_println(ctx);

    /* e^pi; its radius depends on whether the processor has fused
       multiply-add, so print the midpoint and a bound for the radius */
    x = gr_heap_init(ctx);
    status |= gr_pi(x, ctx);
    status |= gr_exp(x, x, ctx);
    gr_ctx_init_real_arb(RR, 200);
    arb_init(y);
    status |= gr_set_other(y, x, ctx, RR);
    arb_printn(y, 40, ARB_STR_NO_RADIUS);
    flint_printf(", radius < 10^-43: %d\n", mag_cmp_2exp_si(arb_radref(y), -143) < 0);
    flint_printf("gamma unsupported: %d\n", gr_gamma(x, x, ctx) == GR_UNABLE);

    gr_mat_init(H, n, n, ctx);
    gr_mat_init(X, n, 1, ctx);
    gr_mat_init(B, n, 1, ctx);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            status |= gr_set_ui(GR_MAT_ENTRY(H, i, j, sz), i + j + 1, ctx);
            status |= gr_inv(GR_MAT_ENTRY(H, i, j, sz), GR_MAT_ENTRY(H, i, j, sz), ctx);
        }
        status |= gr_set_ui(GR_MAT_ENTRY(B, i, 0, sz), 1, ctx);
    }
    status |= gr_mat_nonsingular_solve(X, H, B, ctx);
    status |= gr_mat_print(X, ctx);
    flint_printf("\nstatus = %d\n", status);

    gr_mat_clear(H, ctx);
    gr_mat_clear(X, ctx);
    gr_mat_clear(B, ctx);
    gr_heap_clear(x, ctx);
    arb_clear(y);
    gr_ctx_clear(RR);
    gr_ctx_clear(ctx);
    return 0;
}
```
```output
Real numbers as balls with 3-term double expansion midpoints (d3b)
23.14069263277926900572908636794854738027, radius < 10^-43: 1
gamma unsupported: 1
[[[-6.000000000000000000000000000000000 +/- 1.01e-34]],
[[210.000000000000000000000000000000000 +/- 1.25e-34]],
[[-1680.0000000000000000000000000000000000 +/- 8.02e-35]],
[[5040.0000000000000000000000000000000000 +/- 3.52e-35]],
[[-6300.0000000000000000000000000000000000 +/- 1.14e-35]],
[[2772.00000000000000000000000000000000000 +/- 4.51e-36]]]
status = 0
```

## Choosing a format {#choosing}

- For speed at up to about 106 bits, use `d2` or `d2b`: their
  operations are a few double operations each, and vector operations and
  dot products use SIMD.
- At 159 and 212 bits (`d3`, `d4`), the expansions are still faster
  than `arb_t` but no longer faster than `nfloat` for addition and
  multiplication; they keep an advantage for vectors, division and
  square roots ([](dfloat/performance)).
- Use the ball formats whenever the result must be rigorous; their
  overhead over the plain formats is 1.2 to 1.7 times for products and 2
  to 4 times for division and square roots.
- Above 212 bits, or for functions other than the elementary ones, use
  `arb_t` (rigorous) or `nfloat` (approximate). See
  [](algorithms/floating-point/overview) for a comparison of all
  floating-point types in FLINT.
