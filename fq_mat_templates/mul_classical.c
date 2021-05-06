/*
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

void
TEMPLATE(T, mat_mul_classical) (TEMPLATE(T, mat_t) C,
                                const TEMPLATE(T, mat_t) A,
                                const TEMPLATE(T, mat_t) B,
                                const TEMPLATE(T, ctx_t) ctx)
{
    slong ar, bc, br;
    slong i, j;
    TEMPLATE(T, struct) * trB;
    TMP_INIT;

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
        TEMPLATE(T, mat_swap_entrywise) (C, T, ctx);
        TEMPLATE(T, mat_clear) (T, ctx);
        return;
    }

    TMP_START;

    trB = (TEMPLATE(T, struct) *) TMP_ALLOC(br*bc*sizeof(TEMPLATE(T, struct)));

    /* shallow transpose so columns of B are vectors */
    for (i = 0; i < br; i++)
    {
       for (j = 0; j < bc; j++)
          trB[j*br + i] = *TEMPLATE(T, mat_entry) (B, i, j);
    }
   
    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            _TEMPLATE(T, vec_dot) (TEMPLATE(T, mat_entry) (C, i, j),
               A->rows[i], trB + j * br, br, ctx);

        }
    }

    TMP_END;
}


#endif
