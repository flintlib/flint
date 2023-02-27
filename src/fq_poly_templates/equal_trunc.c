/*
    Copyright (C) 2008, 2009, 2014 William Hart
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
TEMPLATE(T, poly_equal_trunc) (const TEMPLATE(T, poly_t) op1,
                         const TEMPLATE(T, poly_t) op2,
                         slong n, const TEMPLATE(T, ctx_t) ctx)
{
    slong i, len1, len2;

    if (op1 == op2)
        return 1;

    if (n < 0)
       n = 0;

    len1 = FLINT_MIN(op1->length, n);
    len2 = FLINT_MIN(op2->length, n);

    if (len1 < len2)
    {
       for (i = len1; i < len2; i++)
       {
          if (!TEMPLATE(T, is_zero)(op2->coeffs + i, ctx))
             return 0;
       }
    } else if (len2 < len1)
    {
       for (i = len2; i < len1; i++)
       {
          if (!TEMPLATE(T, is_zero)(op1->coeffs + i, ctx))
             return 0;
       }
    }

    for (i = 0; i < FLINT_MIN(len1, len2); i++)
        if (!TEMPLATE(T, equal) (op1->coeffs + i, op2->coeffs + i, ctx))
            return 0;

    return 1;
}


#endif
