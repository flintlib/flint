/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "acb_modular.h"
#include "gr_special.h"
#include "gr_poly.h"
#include "gr_generic.h"

int
gr_generic_hilbert_class_poly(gr_ptr res, slong D, gr_srcptr x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    fmpz_poly_t t;
    fmpz_poly_init(t);

    acb_modular_hilbert_class_poly(t, D);

    if (t->length == 0)
    {
        status = GR_DOMAIN;
    }
    else
    {
        if (ctx->which_ring == GR_CTX_GR_POLY && gr_poly_is_gen(x, POLYNOMIAL_ELEM_CTX(ctx)) == T_TRUE)
            status |= gr_poly_set_fmpz_poly(res, t, POLYNOMIAL_ELEM_CTX(ctx));
        else
            status |= gr_fmpz_poly_evaluate(res, t, x, ctx);
    }

    fmpz_poly_clear(t);

    return status;
}
