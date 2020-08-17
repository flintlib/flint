/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

slong
_TEMPLATE(T, poly_hamming_weight) (const TEMPLATE(T, struct) * op, slong len,
                                   const TEMPLATE(T, ctx_t) ctx)
{
    slong i, sum = 0;
    for (i = 0; i < len; i++)
        sum += !TEMPLATE(T, is_zero) (op + i, ctx);

    return sum;
}

slong
TEMPLATE(T, poly_hamming_weight) (const TEMPLATE(T, poly_t) op,
                                  const TEMPLATE(T, ctx_t) ctx)
{

    return _TEMPLATE(T, poly_hamming_weight) (op->coeffs, op->length, ctx);
}


#endif
