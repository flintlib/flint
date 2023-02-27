/*
    Copyright (C) 2008, 2009 William Hart
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

int
TEMPLATE(T, poly_equal) (const TEMPLATE(T, poly_t) op1,
                         const TEMPLATE(T, poly_t) op2,
                         const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    if (op1 == op2)
        return 1;

    if (op1->length != op2->length)
        return 0;

    for (i = 0; i < op1->length; i++)
        if (!TEMPLATE(T, equal) (op1->coeffs + i, op2->coeffs + i, ctx))
            return 0;

    return 1;
}


#endif
