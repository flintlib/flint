/*
    Copyright (C) 2009, 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fft.h"

void fermat_to_mpz(mpz_t m, mp_limb_t * i, mp_size_t limbs)
{
    mp_limb_signed_t hi;
    mp_ptr mptr = FLINT_MPZ_REALLOC(m, limbs + 1);

    flint_mpn_copyi(mptr, i, limbs + 1);
    hi = i[limbs];
    if (hi < WORD(0))
    {
        mpn_neg(mptr, mptr, limbs + 1);
        m->_mp_size = limbs + 1;
        while ((m->_mp_size) && (!mptr[m->_mp_size - 1]))
            m->_mp_size--;
        m->_mp_size = -m->_mp_size;
    }
    else
    {
        m->_mp_size = limbs + 1;
        while ((m->_mp_size) && (!mptr[m->_mp_size - 1]))
            m->_mp_size--;
    }
}
