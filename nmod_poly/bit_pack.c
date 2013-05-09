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

    Copyright (C) 2007 David Howden
    Copyright (C) 2010 William Hart

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "fmpz.h"

/* Assumes length > 0, bits > 0. */
void
_nmod_poly_bit_pack(mp_ptr res, mp_srcptr poly, len_t len, mp_bitcnt_t bits)
{
    len_t i;
    ulong current_bit = 0, current_limb = 0;
    ulong total_limbs = (len * bits - 1) / FLINT_BITS + 1;
    mp_limb_t temp_lower, temp_upper;

    res[0] = 0L;

    if (bits < FLINT_BITS)
    {
        ulong boundary_limit_bit = FLINT_BITS - bits;

        for (i = 0; i < len; i++)
        {
            if (current_bit > boundary_limit_bit)
            {
                /* the coefficient will be added accross a limb boundary */
                temp_lower = (poly[i] << current_bit);
                temp_upper = (poly[i] >> (FLINT_BITS - current_bit));

                res[current_limb] |= temp_lower;

                current_limb++;
                res[current_limb] = temp_upper;

                current_bit += bits - FLINT_BITS;
            }
            else
            {
                /* the coefficient will fit in the current limb */
                temp_lower = poly[i] << current_bit;
                res[current_limb] |= temp_lower;

                current_bit += bits;

                if (current_bit == FLINT_BITS)
                {
                    current_limb++;
                    if (current_limb < total_limbs)
                        res[current_limb] = 0L;
                    current_bit = 0;
                }
            }
        }
    }
    else if (bits == FLINT_BITS)
    {
        for (i = 0; i < len; i++)
            res[i] = poly[i];
    }
    else if (bits == 2 * FLINT_BITS)
    {
        for (i = 0; i < len; i++)
        {
            res[current_limb++] = poly[i];
            res[current_limb++] = 0L;
        }
    }
    else if (bits < 2 * FLINT_BITS)
    {
        for (i = 0; i < len; i++)
        {
            /* the coefficient will be added accross a limb boundary */
            temp_lower = poly[i] << current_bit;
            temp_upper = r_shift(poly[i], FLINT_BITS - current_bit);

            res[current_limb++] |= temp_lower;
            res[current_limb] = temp_upper;

            current_bit += bits - FLINT_BITS;

            if (current_bit >= FLINT_BITS)
            {
                current_bit -= FLINT_BITS;
                current_limb++;
                if (current_limb < total_limbs)
                    res[current_limb] = 0L;
            }
        }
    }
    else                        /* 2*FLINT_BITS < bits < 3*FLINT_BITS */
    {
        for (i = 0; i < len; i++)
        {
            temp_lower = poly[i] << current_bit;
            temp_upper = r_shift(poly[i], FLINT_BITS - current_bit);

            res[current_limb++] |= temp_lower;
            res[current_limb++] = temp_upper;

            if (current_limb < total_limbs)
                res[current_limb] = 0L;
            current_bit += bits - 2 * FLINT_BITS;

            if (current_bit >= FLINT_BITS)
            {
                current_bit -= FLINT_BITS;
                current_limb++;
                if (current_limb < total_limbs)
                    res[current_limb] = 0L;
            }
        }
    }
}

void
nmod_poly_bit_pack(fmpz_t f, const nmod_poly_t poly,
                   mp_bitcnt_t bit_size)
{
    len_t len, limbs;
    __mpz_struct * mpz;
    len_t i;

    len = nmod_poly_length(poly);

    if (len == 0 || bit_size == 0)
    {
        fmpz_zero(f);
        return;
    }

    mpz = _fmpz_promote(f);
    mpz_realloc2(mpz, len * bit_size);

    limbs = (len * bit_size - 1) / FLINT_BITS + 1;

    _nmod_poly_bit_pack(mpz->_mp_d, poly->coeffs, len, bit_size);

    for (i = limbs - 1; i >= 0; i--)
    {
        if (mpz->_mp_d[i] != 0)
            break;
    }

    mpz->_mp_size = i + 1;
    _fmpz_demote_val(f);
}
