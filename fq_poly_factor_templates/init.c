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
TEMPLATE(T, poly_factor_init) (TEMPLATE(T, poly_factor_t) fac,
                               const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    fac->alloc = 5;
    fac->num = 0;
    fac->poly = flint_malloc(sizeof(TEMPLATE(T, poly_struct)) * 5);
    fac->exp = flint_malloc(sizeof(slong) * 5);

    for (i = 0; i < fac->alloc; i++)
        TEMPLATE(T, poly_init) (fac->poly + i, ctx);
}


#endif
