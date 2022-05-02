/*
    Copyright (C) 2008-2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"

void _fmpz_vec_set_fft(fmpz * coeffs_m, slong length,
                          const mp_ptr * coeffs_f, slong limbs, slong sign)
{
    slong i, size;
    mp_limb_t * data;
    __mpz_struct * mcoeffs_m;

    if (sign)
    {
        for (i = 0; i < length; i++)
        {
            mcoeffs_m = _fmpz_promote(coeffs_m);
            data = FLINT_MPZ_REALLOC(mcoeffs_m, limbs);

			if ((coeffs_f[i][limbs - 1] >> (FLINT_BITS - 1)) || coeffs_f[i][limbs])
            {
                mpn_neg(data, coeffs_f[i], limbs);
                mpn_add_1(data, data, limbs, WORD(1));
                size = limbs;
                while ((size) && (data[size - 1] == 0)) size--; /* normalise */
                mcoeffs_m->_mp_size = -size;
                if (size >= WORD(-1)) _fmpz_demote_val(coeffs_m); /* coefficient may be small*/
            }
            else
            {
                FLINT_MPN_COPYI(data, coeffs_f[i], limbs);
                size = limbs;
                while ((size) && (data[size - 1] == WORD(0))) size--; /* normalise */
                mcoeffs_m->_mp_size = size;
                if (size <= 1) _fmpz_demote_val(coeffs_m); /* coefficient may be small */
            }

            coeffs_m++;
        }
    }
    else
    {
        for (i = 0; i < length; i++)
        {
            mcoeffs_m = _fmpz_promote(coeffs_m);
            data = FLINT_MPZ_REALLOC(mcoeffs_m, limbs);
            FLINT_MPN_COPYI(data, coeffs_f[i], limbs); 
            size = limbs;
            while ((size) && (data[size - 1] == WORD(0))) size--; /* normalise */
            mcoeffs_m->_mp_size = size;
            if (size <= 1) _fmpz_demote_val(coeffs_m); /* coefficient may be small */

            coeffs_m++;
        }
    }
}
