/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void TEMPLATE(T, mat_mul_vec)(TEMPLATE(T, struct) * c,
                                    const TEMPLATE(T, mat_t) A,
                                    const TEMPLATE(T, struct) * b, slong blen,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j;
    slong len = FLINT_MIN(A->c, blen);
    TEMPLATE(T, t) t;

    TEMPLATE(T, init)(t, ctx);

    for (i = A->r - 1; i >= 0; i--)
    {
        TEMPLATE(T, zero)(c + i, ctx);
        for (j = 0; j < len; j++)
        {
            TEMPLATE(T, mul)(t, TEMPLATE(T, mat_entry)(A, i, j), b + j, ctx);
            TEMPLATE(T, add)(c + i, c + i, t, ctx);
        }
    }

    TEMPLATE(T, clear)(t, ctx);
}

void TEMPLATE(T, mat_mul_vec_ptr)(TEMPLATE(T, struct) * const * c,
                            const TEMPLATE(T, mat_t) A,
                            const TEMPLATE(T, struct) * const * b, slong blen,
                            const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j;
    slong len = FLINT_MIN(A->c, blen);
    TEMPLATE(T, t) t;

    TEMPLATE(T, init)(t, ctx);

    for (i = A->r - 1; i >= 0; i--)
    {
        TEMPLATE(T, zero)(c[i], ctx);
        for (j = 0; j < len; j++)
        {
            TEMPLATE(T, mul)(t, TEMPLATE(T, mat_entry)(A, i, j), b[j], ctx);
            TEMPLATE(T, add)(c[i], c[i], t, ctx);
        }
    }

    TEMPLATE(T, clear)(t, ctx);
}

#endif
