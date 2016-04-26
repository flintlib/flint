/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
TEMPLATE(T, mat_mul_classical) (TEMPLATE(T, mat_t) C,
                                const TEMPLATE(T, mat_t) A,
                                const TEMPLATE(T, mat_t) B,
                                const TEMPLATE(T, ctx_t) ctx)
{
    slong ar, bc, br;
    slong i, j, k;
    TEMPLATE(T, t) t;

    ar = A->r;
    br = B->r;
    bc = B->c;

    if (br == 0)
    {
        TEMPLATE(T, mat_zero) (C, ctx);
        return;
    }

    if (C == A || C == B)
    {
        TEMPLATE(T, mat_t) T;
        TEMPLATE(T, mat_init) (T, ar, bc, ctx);
        TEMPLATE(T, mat_mul_classical) (T, A, B, ctx);
        TEMPLATE(T, mat_swap) (C, T, ctx);
        TEMPLATE(T, mat_clear) (T, ctx);
        return;
    }

    TEMPLATE(T, init) (t, ctx);

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            TEMPLATE(T, mul) (TEMPLATE(T, mat_entry) (C, i, j),
                              TEMPLATE(T, mat_entry) (A, i, 0),
                              TEMPLATE(T, mat_entry) (B, 0, j), ctx);

            for (k = 1; k < br; k++)
            {
                TEMPLATE(T, mul) (t,
                                  TEMPLATE(T, mat_entry) (A, i, k),
                                  TEMPLATE(T, mat_entry) (B, k, j), ctx);

                TEMPLATE(T, add) (TEMPLATE(T, mat_entry) (C, i, j),
                                  TEMPLATE(T, mat_entry) (C, i, j), t, ctx);
            }
        }
    }

    TEMPLATE(T, clear) (t, ctx);
}


#endif
