/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

slong
_fmpz_vec_height_index(const fmpz * vec, slong len)
{
    if (len == 1)
    {
        return 0;
    }
    else
    {
        fmpz c;
        mp_srcptr max_d;
        slong max_mpz_limbs, i, max_i, max_coeff, mpz_limbs;

        max_coeff = 0;
        max_i = 0;

        for (i = 0; i < len; i++)
        {
            c = vec[i];

            if (!COEFF_IS_MPZ(c))
            {
                c = FLINT_ABS(c);
                if (c > max_coeff)
                {
                    max_coeff = c;
                    max_i = i;
                }
            }
            else
            {
                __mpz_struct * mc = COEFF_TO_PTR(c);
                max_d = mc->_mp_d;
                max_mpz_limbs = mc->_mp_size;
                max_mpz_limbs = FLINT_ABS(max_mpz_limbs);
                max_i = i;
                i++;
                break;
            }
        }

        for ( ; i < len; i++)
        {
            c = vec[i];

            /* we have found at least one mpz, so only look for those */
            if (COEFF_IS_MPZ(c))
            {
                __mpz_struct * mc = COEFF_TO_PTR(c);
                mpz_limbs = mc->_mp_size;
                mpz_limbs = FLINT_ABS(mpz_limbs);
                if (mpz_limbs > max_mpz_limbs ||
                    ((mpz_limbs == max_mpz_limbs) &&
                    (mpn_cmp(mc->_mp_d, max_d, max_mpz_limbs) > 0)))
                {
                    max_d = mc->_mp_d;
                    max_mpz_limbs = mpz_limbs;
                    max_i = i;
                }
            }
        }

        return max_i;
    }
}
