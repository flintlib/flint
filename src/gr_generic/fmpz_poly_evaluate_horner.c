/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "gr.h"
#include "gr_generic.h"

int
_gr_fmpz_poly_evaluate_horner(gr_ptr res, const fmpz * f, slong len, gr_srcptr x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (len == 0)
    {
        return gr_zero(res, ctx);
    }
    else if (len == 1 || gr_is_zero(x, ctx) == T_TRUE)
    {
        return gr_set_fmpz(res, f, ctx);
    }
    else if (len == 2)
    {
        status |= gr_mul_fmpz(res, x, f + 1, ctx);
        status |= gr_add_fmpz(res, res, f + 0, ctx);
        return status;
    }
    else
    {
        slong i = len - 1;
        gr_ptr t, u;

        GR_TMP_INIT2(t, u, ctx);

        status |= gr_set_fmpz(u, f + i, ctx);

        for (i = len - 2; i >= 0; i--)
        {
            status |= gr_mul(t, u, x, ctx);
            status |= gr_add_fmpz(u, t, f + i, ctx);
        }

        gr_swap(res, u, ctx);

        GR_TMP_CLEAR2(t, u, ctx);

        return status;
    }
}

int
gr_fmpz_poly_evaluate_horner(gr_ptr res, const fmpz_poly_t f, gr_srcptr x, gr_ctx_t ctx)
{
    return _gr_fmpz_poly_evaluate_horner(res, f->coeffs, f->length, x, ctx);
}
