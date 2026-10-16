/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 William Hart
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"
#include "gr_poly.h"

void
_TEMPLATE4(T, poly_evaluate, T, vec_iter)(TEMPLATE(T, struct) * ys,
                                          const TEMPLATE(T, struct) * coeffs, slong len,
                                          const TEMPLATE(T, struct) * xs, slong n,
                                          const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_poly_evaluate_vec_iter(ys, coeffs, len, xs, n, gr_ctx));
}

void
TEMPLATE4(T, poly_evaluate, T, vec_iter)(TEMPLATE(T, struct) * ys,
                                         const TEMPLATE(T, poly_t) poly,
                                         const TEMPLATE(T, struct) * xs, slong n,
                                         const TEMPLATE(T, ctx_t) ctx)
{
    _TEMPLATE4(T, poly_evaluate, T, vec_iter)(ys, poly->coeffs, poly->length,
                                              xs, n, ctx);
}


#endif
