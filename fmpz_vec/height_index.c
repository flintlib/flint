/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

len_t
_fmpz_vec_height_index(const fmpz * vec, len_t len)
{
    if (len == 1)
    {
        return 0;
    }
    else
    {
        fmpz c;
        mp_srcptr max_d;
        len_t max_mpz_limbs, i, max_i, max_coeff, mpz_limbs;

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
                __mpz_struct * mpz_ptr = COEFF_TO_PTR(c);
                max_d = mpz_ptr->_mp_d;
                max_mpz_limbs = mpz_ptr->_mp_size;
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
                __mpz_struct * mpz_ptr = COEFF_TO_PTR(c);
                mpz_limbs = mpz_ptr->_mp_size;
                mpz_limbs = FLINT_ABS(mpz_limbs);
                if (mpz_limbs > max_mpz_limbs ||
                    ((mpz_limbs == max_mpz_limbs) &&
                    (mpn_cmp(mpz_ptr->_mp_d, max_d, max_mpz_limbs) > 0)))
                {
                    max_d = mpz_ptr->_mp_d;
                    max_mpz_limbs = mpz_limbs;
                    max_i = i;
                }
            }
        }

        return max_i;
    }
}
