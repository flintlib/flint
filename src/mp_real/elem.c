/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "impl.h"

void
_mp_real_elem_half_ratio(nn_ptr w, nn_srcptr s, slong n, int minus)
{
    nn_ptr D;
    TMP_INIT;

    if (n <= 2)
    {
        /* the exact floor, in registers */
        ulong s1 = (n == 2) ? s[1] : s[0], s0 = (n == 2) ? s[0] : 0;
        ulong d2, d1, d0, q[2];

        if (minus)
        {
            sub_dddmmmsss(d2, d1, d0, UWORD(2), UWORD(0), UWORD(0), UWORD(0), s1, s0);
        }
        else
        {
            d2 = 2; d1 = s1; d0 = s0;
        }

        _mp_real_divq_4_3z(q, s1, s0, d2, d1, d0);

        if (n == 2)
        {
            w[0] = q[0];
            w[1] = q[1];
        }
        else
            w[0] = q[1];
        return;
    }

    TMP_START;
    D = TMP_ALLOC((n + 1) * sizeof(ulong));
    if (minus)
    {
        /* 2 B^n - s = B^n + (B^n - s) */
        mpn_neg(D, s, n);
        D[n] = 1;
    }
    else
    {
        flint_mpn_copyi(D, s, n);
        D[n] = 2;
    }
    /* floor(s B^n / D) or one more: n limbs, below B^n since s < D */
    flint_mpn_divapprox_fraction(w, s, n, D, n + 1, n);
    TMP_END;
}
