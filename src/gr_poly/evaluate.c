/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

int
_gr_poly_evaluate(gr_ptr y, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t ctx)
{
    return _gr_poly_evaluate_horner(y, f, len, x, ctx);
}

int
gr_poly_evaluate(gr_ptr res, const gr_poly_t f, gr_srcptr a, gr_ctx_t ctx)
{
    return _gr_poly_evaluate(res, f->coeffs, f->length, a, ctx);
}
