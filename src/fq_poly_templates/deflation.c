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

ulong
TEMPLATE(T, poly_deflation) (const TEMPLATE(T, poly_t) input,
                             const TEMPLATE(T, ctx_t) ctx)
{
    ulong deflation;
    slong i, coeff;

    if (input->length <= 1)
        return input->length;

    coeff = 1;
    while (TEMPLATE(T, is_zero) (input->coeffs + coeff, ctx))
        coeff++;

    deflation = n_gcd(input->length - 1, coeff);

    while ((deflation > 1) && (coeff + deflation < input->length))
    {
        for (i = 0; i < deflation - 1; i++)
        {
            coeff++;
            if (!TEMPLATE(T, is_zero) (input->coeffs + coeff, ctx))
                deflation = n_gcd(coeff, deflation);
        }
        if (i == deflation - 1)
            coeff++;
    }

    return deflation;
}


#endif
