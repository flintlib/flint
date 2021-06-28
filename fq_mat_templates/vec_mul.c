/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void TEMPLATE(T, mat_vec_mul)(TEMPLATE(T, struct) * c,
                                    const TEMPLATE(T, struct) * a, slong alen,
                                    const TEMPLATE(T, mat_t) B,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j;
    slong len = FLINT_MIN(B->r, alen);
    TEMPLATE(T, t) t;

    TEMPLATE(T, init)(t, ctx);

    for (i = B->c - 1; i >= 0; i--)
    {
        TEMPLATE(T, zero)(c + i, ctx);
        for (j = 0; j < len; j++)
        {
            TEMPLATE(T, mul)(t, a + j, TEMPLATE(T, mat_entry)(B, j, i), ctx);
            TEMPLATE(T, add)(c + i, c + i, t, ctx);
        }
    }

    TEMPLATE(T, clear)(t, ctx);
}

void TEMPLATE(T, mat_vec_mul_ptr)(TEMPLATE(T, struct) * const * c,
                            const TEMPLATE(T, struct) * const * a, slong alen,
                            const TEMPLATE(T, mat_t) B,
                            const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j;
    slong len = FLINT_MIN(B->r, alen);
    TEMPLATE(T, t) t;

    TEMPLATE(T, init)(t, ctx);

    for (i = B->c - 1; i >= 0; i--)
    {
        TEMPLATE(T, zero)(c[i], ctx);
        for (j = 0; j < len; j++)
        {
            TEMPLATE(T, mul)(t, a[j], TEMPLATE(T, mat_entry)(B, j, i), ctx);
            TEMPLATE(T, add)(c[i], c[i], t, ctx);
        }
    }

    TEMPLATE(T, clear)(t, ctx);
}

#endif
