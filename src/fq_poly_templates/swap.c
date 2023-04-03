/*
    Copyright (C) 2008, 2009 William Hart
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

void
TEMPLATE(T, poly_swap) (TEMPLATE(T, poly_t) op1, TEMPLATE(T, poly_t) op2,
                        const TEMPLATE(T, ctx_t) ctx)
{
    if (op1 != op2)
    {
        slong temp;
        TEMPLATE(T, struct) * temp_c;

        temp = op1->length;
        op1->length = op2->length;
        op2->length = temp;

        temp = op1->alloc;
        op1->alloc = op2->alloc;
        op2->alloc = temp;

        temp_c = op1->coeffs;
        op1->coeffs = op2->coeffs;
        op2->coeffs = temp_c;
    }
}


#endif
