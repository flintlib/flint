/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson
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
TEMPLATE(T, poly_deflate) (TEMPLATE(T, poly_t) result,
                           const TEMPLATE(T, poly_t) input, ulong deflation,
                           const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);

    if (deflation == 0)
    {
        flint_throw(FLINT_DIVZERO, "(%s): Division by zero\n", __func__);
    }

    GR_MUST_SUCCEED(gr_poly_deflate((gr_poly_struct *) result,
                        (const gr_poly_struct *) input, deflation, gr_ctx));
}


#endif
