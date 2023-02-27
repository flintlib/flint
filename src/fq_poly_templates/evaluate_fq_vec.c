/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 William Hart
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
_TEMPLATE4(T, poly_evaluate, T, vec)(TEMPLATE(T, struct) * ys,
                                     const TEMPLATE(T, struct) *coeffs, slong len,
                                     const TEMPLATE(T, struct) *xs, slong n,
                                     const TEMPLATE(T, ctx_t) ctx)
{
    if (len < 32)
        _TEMPLATE4(T, poly_evaluate, T, vec_iter)(ys, coeffs, len,
                                                  xs, n, ctx);
    else
        _TEMPLATE4(T, poly_evaluate, T, vec_fast)(ys, coeffs, len,
                                                  xs, n, ctx);
}

void
TEMPLATE4(T, poly_evaluate, T, vec)(TEMPLATE(T, struct) * ys,
                                    const TEMPLATE(T, poly_t) poly,
                                    const TEMPLATE(T, struct) * xs, slong n,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    _TEMPLATE4(T, poly_evaluate, T, vec)(ys, poly->coeffs, poly->length,
                                         xs, n, ctx);
}


#endif
