/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include "ulong_extras.h"

void
TEMPLATE(T, poly_factor_equal_deg) (TEMPLATE(T, poly_factor_t) factors,
                                    const TEMPLATE(T, poly_t) pol, slong d,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    if (pol->length == d + 1)
    {
        TEMPLATE(T, poly_factor_insert) (factors, pol, 1, ctx);
    }
    else
    {
        TEMPLATE(T, poly_t) f, g, r;
        flint_rand_t state;

        TEMPLATE(T, poly_init) (f, ctx);

        flint_randinit(state);

        while (!TEMPLATE(T, poly_factor_equal_deg_prob)
               (f, state, pol, d, ctx))
        {
        };

        flint_randclear(state);

        TEMPLATE(T, poly_init) (g, ctx);
        TEMPLATE(T, poly_init) (r, ctx);
        TEMPLATE(T, poly_divrem) (g, r, pol, f, ctx);
        TEMPLATE(T, poly_clear) (r, ctx);

        TEMPLATE(T, poly_factor_equal_deg) (factors, f, d, ctx);
        TEMPLATE(T, poly_clear) (f, ctx);
        TEMPLATE(T, poly_factor_equal_deg) (factors, g, d, ctx);
        TEMPLATE(T, poly_clear) (g, ctx);
    }
}


#endif
