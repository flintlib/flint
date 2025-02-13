/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

int
gr_poly_set_str(gr_poly_t res, const char * s, const char * x, gr_ctx_t ctx)
{
    gr_ctx_t rctx;
    int status = GR_SUCCESS;
    gr_ctx_init_gr_poly(rctx, ctx);
    status |= gr_ctx_set_gen_name(rctx, x);
    status |= gr_set_str(res, s, rctx);
    gr_ctx_clear(rctx);
    return status;
}

int
_gr_poly_set_str(gr_ptr res, const char * s, const char * x, slong len, gr_ctx_t ctx)
{
    gr_poly_t t;
    int status;
    gr_poly_init(t, ctx);
    status = gr_poly_set_str(t, s, x, ctx);

    if (status == GR_SUCCESS)
    {
        if (t->length <= len)
        {
            _gr_vec_swap(res, t->coeffs, t->length, ctx);
            status |= _gr_vec_zero(GR_ENTRY(res, t->length, ctx->sizeof_elem), len - t->length, ctx);
            /* make sure t is a valid polynomial before clearing
               just so that we don't strike the debug assertion in gr_poly_clear */
            _gr_poly_set_length(t, 0, ctx);
        }
        else
        {
            status |= _gr_vec_zero(res, len, ctx);
            status = GR_UNABLE;
        }
    }

    gr_poly_clear(t, ctx);
    return status;
}
