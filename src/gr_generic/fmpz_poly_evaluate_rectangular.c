/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz_poly.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_vec.h"

int
_gr_fmpz_poly_evaluate_rectangular(gr_ptr res, const fmpz * poly, slong len, gr_srcptr x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i, m, r;
    gr_ptr xs, s, t, c;
    slong sz = ctx->sizeof_elem;

    if (len <= 2)
        return _gr_fmpz_poly_evaluate_horner(res, poly, len, x, ctx);

    m = n_sqrt(len) + 1;
    r = (len + m - 1) / m;

    GR_TMP_INIT_VEC(xs, m + 1, ctx);
    GR_TMP_INIT3(s, t, c, ctx);

    status = _gr_vec_set_powers(xs, x, m + 1, ctx);

    status |= gr_set_fmpz(res, poly + (r - 1) * m, ctx);
    status |= _gr_vec_dot_fmpz(res, res, 0, GR_ENTRY(xs, 1, sz), poly + (r - 1) * m + 1, len - (r - 1) * m - 1, ctx);

    for (i = r - 2; i >= 0; i--)
    {
        status |= gr_set_fmpz(s, poly + i * m, ctx);
        status |= _gr_vec_dot_fmpz(s, s, 0, GR_ENTRY(xs, 1, sz), poly + i * m + 1, m - 1, ctx);
        status |= gr_mul(res, res, GR_ENTRY(xs, m, sz), ctx);
        status |= gr_add(res, res, s, ctx);
    }

    GR_TMP_CLEAR_VEC(xs, m + 1, ctx);
    GR_TMP_CLEAR3(s, t, c, ctx);

    return status;
}

int
gr_fmpz_poly_evaluate_rectangular(gr_ptr res, const fmpz_poly_t f, gr_srcptr x, gr_ctx_t ctx)
{
    return _gr_fmpz_poly_evaluate_rectangular(res, f->coeffs, f->length, x, ctx);
}
