/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_extras.h"

void
fmpz_lshift_mpn(fmpz_t z, mp_srcptr d, mp_size_t dn, int sgnbit, flint_bitcnt_t shift)
{
    __mpz_struct * zmpz;
    mp_ptr zp;
    mp_size_t zn, shift_limbs;
    flint_bitcnt_t shift_bits;

    zmpz = _fmpz_promote(z);

    shift_limbs = shift / FLINT_BITS;
    shift_bits = shift % FLINT_BITS;
    zn = dn + shift_limbs + (shift_bits != 0);

    if (zmpz->_mp_alloc < zn)
        mpz_realloc2(zmpz, zn * FLINT_BITS);

    zp = zmpz->_mp_d;
    flint_mpn_zero(zp, shift_limbs);

    if (shift_bits == 0)
    {
        flint_mpn_copyi(zp + shift_limbs, d, dn);
    }
    else
    {
        zp[zn - 1] = mpn_lshift(zp + shift_limbs, d, dn, shift_bits);
        while (zp[zn - 1] == 0)
            zn--;
    }

    zmpz->_mp_size = sgnbit ? -(slong) zn : zn;
    _fmpz_demote_val(z);
}

