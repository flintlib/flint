/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2012 Andres Goens
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
TEMPLATE(T, poly_randtest) (TEMPLATE(T, poly_t) f, flint_rand_t state,
                            slong len, const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    TEMPLATE(T, poly_fit_length) (f, len, ctx);
    for (i = 0; i < len; i++)
    {
        TEMPLATE(T, randtest) (f->coeffs + i, state, ctx);
    }
    _TEMPLATE(T, poly_set_length) (f, len, ctx);
    _TEMPLATE(T, poly_normalise) (f, ctx);
}

void
TEMPLATE(T, poly_randtest_not_zero) (TEMPLATE(T, poly_t) f, flint_rand_t state,
                                     slong len, const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    if (len == 0)
    {
        flint_throw(FLINT_ERROR, "(%s): len = 0\n", __func__);
    }

    TEMPLATE(T, poly_randtest) (f, state, len, ctx);
    for (i = 0; (i < 10) && TEMPLATE(T, poly_is_zero) (f, ctx); i++)
        TEMPLATE(T, poly_randtest) (f, state, len, ctx);
    if (TEMPLATE(T, poly_is_zero) (f, ctx))
        TEMPLATE(T, poly_one) (f, ctx);
}


#endif
