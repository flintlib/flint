/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_dft.h"

/* Sets res = w^e * x (or w^(-e) * x if inverse is set), where w is the
   root of unity of the plan P and 0 <= e < n. res must not alias x.

   In complex mode, rtmp must point to an initialized temporary element
   of P->real_ctx. Roots at the quarter and eighth points, known
   structurally since complex mode uses the standard root
   w = exp(-2 pi i / n), are handled with 0 or 2 real multiplications.
   A general root c + d*i is handled with a single multiplication in
   the complex ring, or, when the Karatsuba tables are enabled, with
   3 real multiplications using the precomputed entries
   (c, d - c, d + c):

       k1 = c * (xr + xi)
       re = k1 - xi * (d + c) = c * xr - d * xi
       im = k1 + xr * (d - c) = c * xi + d * xr
*/
int
_gr_dft_mul_root(gr_ptr res, gr_srcptr x, ulong e, int inverse,
        gr_ptr rtmp, const gr_dft_pre_t P)
{
    int status = GR_SUCCESS;
    gr_ctx_struct * rctx;
    slong rsz;
    gr_srcptr xr, xi, c;
    gr_ptr rr, ri;

    if (inverse && e != 0)
        e = P->n - e;

    if (e == 0)
        return gr_set(res, x, P->ctx);

    if (P->real_ctx == NULL)
        return gr_mul(res, x,
                GR_ENTRY(P->roots, e, P->ctx->sizeof_elem), P->ctx);

    rctx = P->real_ctx;
    rsz = rctx->sizeof_elem;
    xr = x;
    xi = GR_ENTRY(x, 1, rsz);
    rr = res;
    ri = GR_ENTRY(res, 1, rsz);
    /* real part of the root, valid whether or not wtab is present */
    c = GR_ENTRY(P->roots, e, P->ctx->sizeof_elem);

    switch (P->wclass[e])
    {
        case GR_DFT_ROOT_NEG_ONE:
            /* res = -x */
            status |= gr_neg(rr, xr, rctx);
            status |= gr_neg(ri, xi, rctx);
            break;

        case GR_DFT_ROOT_I:
            /* i x = -xi + xr i */
            status |= gr_neg(rr, xi, rctx);
            status |= gr_set(ri, xr, rctx);
            break;

        case GR_DFT_ROOT_NEG_I:
            /* -i x = xi - xr i */
            status |= gr_set(rr, xi, rctx);
            status |= gr_neg(ri, xr, rctx);
            break;

        case GR_DFT_ROOT_DMC_ZERO:
            /* c (1 + i) x = c (xr - xi) + c (xr + xi) i */
            status |= gr_sub(rtmp, xr, xi, rctx);
            status |= gr_mul(rr, c, rtmp, rctx);
            status |= gr_add(rtmp, xr, xi, rctx);
            status |= gr_mul(ri, c, rtmp, rctx);
            break;

        case GR_DFT_ROOT_DPC_ZERO:
            /* c (1 - i) x = c (xr + xi) + c (xi - xr) i */
            status |= gr_add(rtmp, xr, xi, rctx);
            status |= gr_mul(rr, c, rtmp, rctx);
            status |= gr_sub(rtmp, xi, xr, rctx);
            status |= gr_mul(ri, c, rtmp, rctx);
            break;

        default:
            if (P->wtab != NULL)
            {
                /* 3-multiplication complex Karatsuba product */
                status |= _gr_dft_mul_const(res, x, NULL,
                        GR_ENTRY(P->wtab, 3 * e, rsz), rtmp, P);
            }
            else
            {
                status |= gr_mul(res, x,
                        GR_ENTRY(P->roots, e, P->ctx->sizeof_elem), P->ctx);
            }
            break;
    }

    return status;
}

/* Sets res = w * x for a fixed plan constant w, given the precomputed
   complex Karatsuba table entry wtab3 = (c, d - c, d + c) if it is not
   NULL (in which case w may be NULL, and rtmp must point to an
   initialized temporary element of P->real_ctx), or the constant w
   itself with wtab3 == NULL otherwise. res must not alias x. */
int
_gr_dft_mul_const(gr_ptr res, gr_srcptr x, gr_srcptr w, gr_srcptr wtab3,
        gr_ptr rtmp, const gr_dft_pre_t P)
{
    int status = GR_SUCCESS;
    gr_ctx_struct * rctx;
    slong rsz;
    gr_srcptr xr, xi, c, dmc, dpc;
    gr_ptr rr, ri;

    if (wtab3 == NULL)
        return gr_mul(res, x, w, P->ctx);

    rctx = P->real_ctx;
    rsz = rctx->sizeof_elem;
    xr = x;
    xi = GR_ENTRY(x, 1, rsz);
    rr = res;
    ri = GR_ENTRY(res, 1, rsz);
    c = wtab3;
    dmc = GR_ENTRY(c, 1, rsz);
    dpc = GR_ENTRY(c, 2, rsz);

    status |= gr_add(rtmp, xr, xi, rctx);
    status |= gr_mul(rtmp, c, rtmp, rctx);      /* k1 */
    status |= gr_mul(rr, xi, dpc, rctx);        /* k3 */
    status |= gr_mul(ri, xr, dmc, rctx);        /* k2 */
    status |= gr_sub(rr, rtmp, rr, rctx);
    status |= gr_add(ri, rtmp, ri, rctx);

    return status;
}
