/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz
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
TEMPLATE(T, poly_factor_fit_length) (TEMPLATE(T, poly_factor_t) fac, slong len,
                                     const TEMPLATE(T, ctx_t) ctx)
{
    if (len > fac->alloc)
    {
        /* At least double number of allocated coeffs */
        if (len < 2 * fac->alloc)
            len = 2 * fac->alloc;
        TEMPLATE(T, poly_factor_realloc) (fac, len, ctx);
    }
}


#endif
