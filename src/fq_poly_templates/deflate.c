/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson
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
TEMPLATE(T, poly_deflate) (TEMPLATE(T, poly_t) result,
                           const TEMPLATE(T, poly_t) input, ulong deflation,
                           const TEMPLATE(T, ctx_t) ctx)
{
    slong res_length, i;

    if (deflation == 0)
    {
        flint_throw(FLINT_ERROR, "(%s): Division by zero\n", __func__);
    }

    if (input->length <= 1 || deflation == 1)
    {
        TEMPLATE(T, poly_set) (result, input, ctx);
        return;
    }

    res_length = (input->length - 1) / deflation + 1;
    TEMPLATE(T, poly_fit_length) (result, res_length, ctx);
    for (i = 0; i < res_length; i++)
        TEMPLATE(T, set) (result->coeffs + i, input->coeffs + (i * deflation),
                          ctx);

    result->length = res_length;
}


#endif
