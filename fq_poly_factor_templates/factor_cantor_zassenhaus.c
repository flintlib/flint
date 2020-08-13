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

void
TEMPLATE(T, poly_factor_cantor_zassenhaus) (TEMPLATE(T, poly_factor_t) res,
                                            const TEMPLATE(T, poly_t) f,
                                            const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_t) h, v, g, x;
    fmpz_t q;
    slong i, j, num;

    fmpz_init(q);
    TEMPLATE(T, ctx_order) (q, ctx);

    TEMPLATE(T, poly_init) (h, ctx);
    TEMPLATE(T, poly_init) (g, ctx);
    TEMPLATE(T, poly_init) (v, ctx);
    TEMPLATE(T, poly_init) (x, ctx);

    TEMPLATE(T, poly_gen) (h, ctx);
    TEMPLATE(T, poly_gen) (x, ctx);

    TEMPLATE(T, poly_make_monic) (v, f, ctx);

    i = 0;
    do
    {
        i++;
        TEMPLATE(T, poly_powmod_fmpz_binexp) (h, h, q, v, ctx);

        TEMPLATE(T, poly_sub) (h, h, x, ctx);
        TEMPLATE(T, poly_gcd) (g, h, v, ctx);
        TEMPLATE(T, poly_add) (h, h, x, ctx);

        if (g->length != 1)
        {
            TEMPLATE(T, poly_make_monic) (g, g, ctx);
            num = res->num;

            TEMPLATE(T, poly_factor_equal_deg) (res, g, i, ctx);
            for (j = num; j < res->num; j++)
                res->exp[j] = TEMPLATE(T, poly_remove) (v, res->poly + j, ctx);
        }
    }
    while (v->length >= 2 * i + 3);

    if (v->length > 1)
        TEMPLATE(T, poly_factor_insert) (res, v, 1, ctx);

    TEMPLATE(T, poly_clear) (g, ctx);
    TEMPLATE(T, poly_clear) (h, ctx);
    TEMPLATE(T, poly_clear) (v, ctx);
    TEMPLATE(T, poly_clear) (x, ctx);
    fmpz_clear(q);
}


#endif
