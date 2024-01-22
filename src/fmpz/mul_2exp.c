/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2022, 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "fmpz.h"

void
_fmpz_mul_2exp(fmpz_t rp, fmpz xs, ulong exp, const ulong * xd)
{
    mpz_ptr rm;
    mpz_srcptr xm;
    slong rsz, xabssz, xsz;
    mp_ptr rd;
    slong limbshift;

    rm = _fmpz_promote(rp);

    if (!COEFF_IS_MPZ(xs))
    {
        xm = NULL;
        xsz = xs;
        xabssz = 1;
    }
    else
    {
        xm = COEFF_TO_PTR(xs);
        xd = xm->_mp_d;
        xsz = xm->_mp_size;
        xabssz = FLINT_ABS(xsz);
    }

    limbshift = exp / FLINT_BITS;
    exp %= FLINT_BITS;

    rsz = xabssz + 1 + limbshift;

    if (rm->_mp_alloc < rsz)
    {
        _mpz_realloc(rm, rsz);
        if (xm == rm)
            xd = xm->_mp_d;
    }

    rd = rm->_mp_d;

    if (exp != 0)
        rd[rsz - 1] = mpn_lshift(rd + limbshift, xd, xabssz, exp);
    else
        memmove(rd + limbshift, xd, sizeof(mp_limb_t) * xabssz);

    for (slong ix = 0; ix < limbshift; ix++)
        rd[ix] = 0;

    rm->_mp_size = rsz - (rd[rsz - 1] == 0 || exp == 0);
    if (xsz < 0)
        rm->_mp_size = -rm->_mp_size;
}
