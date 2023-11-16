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


/* todo: divconquer algorithm */
/* todo: algorithm selection */
int
_gr_fmpz_poly_evaluate(gr_ptr res, const fmpz * f, slong len, gr_srcptr x, gr_ctx_t ctx)
{
    if (len <= 6)
        return _gr_fmpz_poly_evaluate_horner(res, f, len, x, ctx);
    else
        return _gr_fmpz_poly_evaluate_rectangular(res, f, len, x, ctx);
}

int
gr_fmpz_poly_evaluate(gr_ptr res, const fmpz_poly_t f, gr_srcptr x, gr_ctx_t ctx)
{
    return _gr_fmpz_poly_evaluate(res, f->coeffs, f->length, x, ctx);
}
