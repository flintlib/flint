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
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
TEMPLATE(T, poly_factor_insert) (TEMPLATE(T, poly_factor_t) fac,
                                 const TEMPLATE(T, poly_t) poly, slong exp,
                                 const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    if (poly->length <= 1)
        return;

    for (i = 0; i < fac->num; i++)
    {
        if (TEMPLATE(T, poly_equal) (poly, fac->poly + i, ctx))
        {
            fac->exp[i] += exp;
            return;
        }
    }

    if (fac->alloc == fac->num)
    {
        slong new_size = 2 * fac->alloc;

        fac->poly =
            flint_realloc(fac->poly,
                          sizeof(TEMPLATE(T, poly_struct)) * new_size);
        fac->exp = flint_realloc(fac->exp, sizeof(slong) * new_size);

        for (i = fac->alloc; i < new_size; i++)
            TEMPLATE(T, poly_init) (fac->poly + i, ctx);

        fac->alloc = new_size;
    }

    TEMPLATE(T, poly_set) (fac->poly + fac->num, poly, ctx);
    fac->exp[fac->num] = exp;
    fac->num++;
}


#endif
