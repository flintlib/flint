/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

void fmpz_get_signed_ui_array(mp_limb_t * r, slong n, const fmpz_t x)
{
    int neg;
    slong i, sz;

    FLINT_ASSERT(n > 0);

    if (!COEFF_IS_MPZ(*x))
    {
        neg = *x < 0;
        r[0] = FLINT_ABS(*x);
        i = 1;
    }
    else
    {
        __mpz_struct * p = COEFF_TO_PTR(*x);
        neg = p->_mp_size < 0;
        sz = FLINT_ABS(p->_mp_size);

        for (i = 0; i < n && i < sz; i++)
            r[i] = p->_mp_d[i];
    }

    for ( ; i < n; i++)
        r[i] = 0;

    if (neg)
        mpn_neg(r, r, n);
}

