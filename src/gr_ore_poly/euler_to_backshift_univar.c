/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Opus 4.8 */

#include "gr.h"
#include "gr_poly.h"
#include "gr_ore_poly.h"

/* Bridge from the Euler operator (K[x]<theta>, theta = x d/dx) to the backward
   shift operator (K[n]<S^{-1}>, S^{-1} a(n) = a(n-1)) under the exact
   isomorphism theta <-> n, x <-> S^{-1}. Writing op = sum_k x^k Q_k(theta), the
   image is sum_k Q_k(n-k) S^{-k}: transpose the (theta-degree, x-degree)
   coefficient matrix and Taylor-shift column k by -k. The caller must allocate
   res to reslen, the number of S^{-1}-coefficients (= 1 + max x-degree of
   op). */
int
_gr_ore_poly_euler_to_backshift_univar(gr_ptr res, slong reslen, gr_srcptr op, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ctx_struct * sctx;
    slong bsz = ctx->sizeof_elem, ssz;
    slong i, k;
    gr_ptr c;

    if (ctx->which_ring != GR_CTX_GR_POLY)
        return GR_UNABLE;

    sctx = POLYNOMIAL_ELEM_CTX(ctx);
    ssz = sctx->sizeof_elem;

    GR_TMP_INIT(c, sctx);
    for (k = 0; k < reslen; k++)
    {
        gr_poly_struct * rk = (gr_poly_struct *) GR_ENTRY(res, k, bsz);

        gr_poly_fit_length(rk, len, sctx);
        for (i = 0; i < len; i++)
        {
            const gr_poly_struct * oi = (const gr_poly_struct *) GR_ENTRY(op, i, bsz);
            gr_ptr dst = GR_ENTRY(rk->coeffs, i, ssz);
            if (k < oi->length)
                status |= gr_set(dst, GR_ENTRY(oi->coeffs, k, ssz), sctx);
            else
                status |= gr_zero(dst, sctx);
        }
        _gr_poly_set_length(rk, len, sctx);
        _gr_poly_normalise(rk, sctx);

        status |= gr_set_si(c, -k, sctx);
        status |= gr_poly_taylor_shift(rk, rk, c, sctx);
    }
    GR_TMP_CLEAR(c, sctx);

    return status;
}
