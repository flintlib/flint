/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arf.h"

int
arf_dot(arf_t res, const arf_t initial, int subtract,
    arf_srcptr x, slong xstep, arf_srcptr y, slong ystep, slong len, slong prec, arf_rnd_t rnd)
{
    arf_ptr A;
    slong i;
    int inexact;

    if (len <= 0)
    {
        if (initial == NULL)
        {
            arf_zero(res);
            return 0;
        }

        return arf_set_round(res, initial, prec, rnd);
    }

    A = _arf_vec_init(len + (initial != NULL));

    for (i = 0; i < len; i++)
    {
        arf_mul(A + i, x + i * xstep, y + i * ystep, ARF_PREC_EXACT, ARF_RND_DOWN);
        if (subtract)
            arf_neg(A + i, A + i);
    }

    if (initial != NULL)
        arf_set(A + len, initial);

    inexact = arf_sum(res, A, len + (initial != NULL), prec, rnd);

    _arf_vec_clear(A, len + (initial != NULL));

    return inexact;
}
