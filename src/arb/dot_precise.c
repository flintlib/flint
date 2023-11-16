/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_dot_precise(arb_t res, const arb_t initial, int subtract,
    arb_srcptr x, slong xstep, arb_srcptr y, slong ystep, slong len, slong prec)
{
    arf_ptr A, B;
    arf_t t, u;
    slong i;
    int inexact;

    if (len <= 0)
    {
        if (initial == NULL)
            arb_zero(res);
        else
            arb_set_round(res, initial, prec);
        return;
    }

    A = flint_calloc(len + (initial != NULL), sizeof(arf_struct));
    B = flint_calloc(3 * len + 1 + (initial != NULL), sizeof(arf_struct));

    for (i = 0; i < len; i++)
    {
        arb_srcptr xp = x + i * xstep;
        arb_srcptr yp = y + i * ystep;

        arf_mul(A + i, arb_midref(xp), arb_midref(yp), ARF_PREC_EXACT, ARF_RND_DOWN);
        if (subtract)
            arf_neg(A + i, A + i);

        arf_init_set_mag_shallow(t, arb_radref(xp));
        arf_init_set_mag_shallow(u, arb_radref(yp));

        arf_mul(B + 3 * i, t, u, ARF_PREC_EXACT, ARF_RND_DOWN);

        arf_mul(B + 3 * i + 1, t, arb_midref(yp), ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_abs(B + 3 * i + 1, B + 3 * i + 1);

        arf_mul(B + 3 * i + 2, u, arb_midref(xp), ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_abs(B + 3 * i + 2, B + 3 * i + 2);
    }

    if (initial != NULL)
    {
        arf_set(A + len, arb_midref(initial));
        arf_set_mag(B + 3 * len + 1, arb_radref(initial));
    }

    inexact = arf_sum(arb_midref(res), A, len + (initial != NULL), prec, ARF_RND_DOWN);
    if (inexact)
        arf_mag_set_ulp(arb_radref(res), arb_midref(res), prec);
    else
        mag_zero(arb_radref(res));
    arf_set_mag(B + 3 * len, arb_radref(res));

    arf_sum(A, B, 3 * len + 1 + (initial != NULL), 3 * MAG_BITS, ARF_RND_UP);
    arf_get_mag(arb_radref(res), A);

    for (i = 0; i < len + (initial != NULL); i++)
        arf_clear(A + i);
    for (i = 0; i < 3 * len + 1 + (initial != NULL); i++)
        arf_clear(B + i);

    flint_free(A);
    flint_free(B);
}
