/*
    Copyright (C) 2026 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
TEMPLATE(T, mat_scalar_mul_fmpz)(TEMPLATE(T, mat_t) B,
                          const TEMPLATE(T, mat_t) A,
                          const fmpz_t x,
                          const TEMPLATE(T, ctx_t) ctx)
{
    slong i;

    if (B->c < 1)
        return;

    for (i = 0; i < B->r; i++)
        _TEMPLATE(T, vec_scalar_mul_fmpz) (
            TEMPLATE(T, mat_entry)(B, i, 0),
            TEMPLATE(T, mat_entry)(A, i, 0),
            B->c,
            x,
            ctx);
}


#endif
