/*
    Copyright (C) 2018 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T
#include "templates.h"

int
TEMPLATE(T, mat_inv)(TEMPLATE(T, mat_t) B, TEMPLATE(T, mat_t) A,
                                   const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, mat_t) I;
    slong i, dim;
    int result;

    dim = A->r;

    switch (dim)
    {
        case 0:
            result = 1;
            break;

        case 1:
            if (TEMPLATE(T, is_zero)(TEMPLATE(T, mat_entry)(A, 0, 0), ctx))
            {
                result = 0;
            }
            else
            {
                TEMPLATE(T, inv)(TEMPLATE(T, mat_entry(B, 0, 0)), TEMPLATE(T, mat_entry(A, 0, 0)), ctx);
                result = 1;
            }
            break;

        default:
            TEMPLATE(T, mat_init)(I, dim, dim, ctx);
            for (i = 0; i < dim; i++)
                TEMPLATE(T, one)(TEMPLATE(T, mat_entry)(I, i, i), ctx);
            result = TEMPLATE(T, mat_solve)(B, A, I, ctx);
            TEMPLATE(T, mat_clear)(I, ctx);
    }

    return result;
}

#endif
