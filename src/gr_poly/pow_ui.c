/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

int
_gr_poly_pow_ui(gr_ptr res,
    gr_srcptr f, slong flen, ulong exp, gr_ctx_t ctx)
{
    return _gr_poly_pow_ui_binexp(res, f, flen, exp, ctx);
}

int
gr_poly_pow_ui(gr_poly_t res,
    const gr_poly_t poly, ulong exp, gr_ctx_t ctx)
{
    return gr_poly_pow_ui_binexp(res, poly, exp, ctx);
}
