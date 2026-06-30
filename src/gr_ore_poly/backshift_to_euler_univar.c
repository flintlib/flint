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

/* Bridge from the backward shift operator (K[n]<S^{-1}>) to the Euler operator
   (K[x]<theta>), the inverse of _gr_ore_poly_euler_to_backshift: Taylor-shift
   column k by +k, then transpose the (S^{-1}-degree, n-degree) coefficient
   matrix. The caller must allocate res to reslen, the number of
   theta-coefficients (= 1 + max n-degree of op). */
int
_gr_ore_poly_backshift_to_euler_univar(gr_ptr res, slong reslen, gr_srcptr op, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ctx_struct * sctx;
    slong bsz = ctx->sizeof_elem, ssz;
    slong i, k;
    gr_ptr c, cols;

    if (ctx->which_ring != GR_CTX_GR_POLY)
        return GR_UNABLE;

    sctx = POLYNOMIAL_ELEM_CTX(ctx);
    ssz = sctx->sizeof_elem;

    GR_TMP_INIT(c, sctx);
    GR_TMP_INIT_VEC(cols, len, ctx);
    for (k = 0; k < len; k++)
    {
        const gr_poly_struct * ok = (const gr_poly_struct *) GR_ENTRY(op, k, bsz);
        gr_poly_struct * colk = (gr_poly_struct *) GR_ENTRY(cols, k, bsz);
        status |= gr_poly_set(colk, ok, sctx);
        status |= gr_set_si(c, k, sctx);
        status |= gr_poly_taylor_shift(colk, colk, c, sctx);
    }

    for (i = 0; i < reslen; i++)
    {
        gr_poly_struct * ri = (gr_poly_struct *) GR_ENTRY(res, i, bsz);

        gr_poly_fit_length(ri, len, sctx);
        for (k = 0; k < len; k++)
        {
            gr_poly_struct * colk = (gr_poly_struct *) GR_ENTRY(cols, k, bsz);
            gr_ptr dst = GR_ENTRY(ri->coeffs, k, ssz);
            if (i < colk->length)
                status |= gr_set(dst, GR_ENTRY(colk->coeffs, i, ssz), sctx);
            else
                status |= gr_zero(dst, sctx);
        }
        _gr_poly_set_length(ri, len, sctx);
        _gr_poly_normalise(ri, sctx);
    }

    GR_TMP_CLEAR_VEC(cols, len, ctx);
    GR_TMP_CLEAR(c, sctx);

    return status;
}
