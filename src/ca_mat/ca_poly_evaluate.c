/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"
#include "gr.h"
#include "gr_mat.h"

void
_ca_mat_ca_poly_evaluate(ca_mat_t y, ca_srcptr poly,
    slong len, const ca_mat_t x, ca_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_ca_from_ref(gr_ctx, GR_CTX_CC_CA, ctx);
    GR_MUST_SUCCEED(_gr_mat_gr_poly_evaluate((gr_mat_struct *) y, poly, len,
        (const gr_mat_struct *) x, gr_ctx));
}

void
ca_mat_ca_poly_evaluate(ca_mat_t res, const ca_poly_t f, const ca_mat_t a, ca_ctx_t ctx)
{
    _ca_mat_ca_poly_evaluate(res, f->coeffs, f->length, a, ctx);
}
