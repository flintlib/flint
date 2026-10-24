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

void TEMPLATE(T, mat_mul_vec)(TEMPLATE(T, struct) * c,
                                    const TEMPLATE(T, mat_t) A,
                                    const TEMPLATE(T, struct) * b, slong blen,
                                    const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    slong len = FLINT_MAX(FLINT_MIN(A->c, blen), 0);

    for (i = 0; i < A->r; i++)
        _TEMPLATE(T, vec_dot)(c + i, TEMPLATE(T, mat_entry)(A, i, 0), b, len, ctx);
}

void TEMPLATE(T, mat_mul_vec_ptr)(TEMPLATE(T, struct) * const * c,
                            const TEMPLATE(T, mat_t) A,
                            const TEMPLATE(T, struct) * const * b, slong blen,
                            const TEMPLATE(T, ctx_t) ctx)
{
    slong i;
    slong len = FLINT_MAX(FLINT_MIN(A->c, blen), 0);
    TEMPLATE(T, struct) * bb;
    TMP_INIT;

    TMP_START;

    bb = TMP_ARRAY_ALLOC(len, TEMPLATE(T, struct));

    /* shallow copies */
    for (i = 0; i < len; i++)
        bb[i] = *b[i];

    for (i = 0; i < A->r; i++)
        _TEMPLATE(T, vec_dot)(c[i], TEMPLATE(T, mat_entry)(A, i, 0), bb, len, ctx);

    TMP_END;
}

#endif
