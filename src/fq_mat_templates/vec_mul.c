/*
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2026 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
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
    slong len = FLINT_MAX(FLINT_MIN(B->r, alen), 0);
    TEMPLATE(T, struct) * col;
    TMP_INIT;

    TMP_START;

    col = TMP_ARRAY_ALLOC(len, TEMPLATE(T, struct));

    for (j = 0; j < B->c; j++)
    {
        /* shallow copies */
        for (i = 0; i < len; i++)
            col[i] = *TEMPLATE(T, mat_entry)(B, i, j);

        _TEMPLATE(T, vec_dot)(c + j, a, col, len, ctx);
    }

    TMP_END;
}

void TEMPLATE(T, mat_vec_mul_ptr)(TEMPLATE(T, struct) * const * c,
                            const TEMPLATE(T, struct) * const * a, slong alen,
                            const TEMPLATE(T, mat_t) B,
                            const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j;
    slong len = FLINT_MAX(FLINT_MIN(B->r, alen), 0);
    TEMPLATE(T, struct) * aa, * col;
    TMP_INIT;

    TMP_START;

    aa = TMP_ARRAY_ALLOC(2 * len, TEMPLATE(T, struct));
    col = aa + len;

    /* shallow copies */
    for (i = 0; i < len; i++)
        aa[i] = *a[i];

    for (j = 0; j < B->c; j++)
    {
        for (i = 0; i < len; i++)
            col[i] = *TEMPLATE(T, mat_entry)(B, i, j);

        _TEMPLATE(T, vec_dot)(c[j], aa, col, len, ctx);
    }

    TMP_END;
}

#endif
