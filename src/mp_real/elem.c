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
