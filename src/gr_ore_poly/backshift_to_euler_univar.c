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

int
_gr_ore_poly_backshift_to_euler_univar(gr_ptr res, slong reslen, gr_srcptr op, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    gr_ctx_struct * sctx;
    slong bsz = ctx->sizeof_elem, ssz;
    slong i, k;
    gr_ptr _k, rcoeffs;

    if (ctx->which_ring != GR_CTX_GR_POLY)
        return GR_UNABLE;

    sctx = POLYNOMIAL_ELEM_CTX(ctx);
    ssz = sctx->sizeof_elem;

    GR_TMP_INIT(_k, sctx);
    GR_TMP_INIT_VEC(rcoeffs, len, ctx);

    for (k = 0; k < len; k++)
    {
        gr_poly_struct * rck = (gr_poly_struct *) GR_ENTRY(rcoeffs, k, bsz);
        status |= gr_poly_set(rck, GR_ENTRY(op, k, bsz), sctx);
        status |= gr_set_si(_k, k, sctx);
        status |= gr_poly_taylor_shift(rck, rck, _k, sctx);
    }

    for (i = 0; i < reslen; i++)
    {
        gr_poly_struct * ri = (gr_poly_struct *) GR_ENTRY(res, i, bsz);

        gr_poly_fit_length(ri, len, sctx);
        for (k = 0; k < len; k++)
        {
            gr_poly_struct * rck = (gr_poly_struct *) GR_ENTRY(rcoeffs, k, bsz);
            gr_ptr dst = GR_ENTRY(ri->coeffs, k, ssz);
            if (i < rck->length)
                status |= gr_set(dst, GR_ENTRY(rck->coeffs, i, ssz), sctx);
            else
                status |= gr_zero(dst, sctx);
        }
        _gr_poly_set_length(ri, len, sctx);
        _gr_poly_normalise(ri, sctx);
    }

    GR_TMP_CLEAR_VEC(rcoeffs, len, ctx);
    GR_TMP_CLEAR(_k, sctx);

    return status;
}
