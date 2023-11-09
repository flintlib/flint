/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpz.h"

int
fmpz_bit_pack(mp_ptr arr, flint_bitcnt_t shift, flint_bitcnt_t bits,
              const fmpz_t coeff, int negate, int borrow)
{
    mp_limb_t save = arr[0];
    fmpz c = *coeff;
    int sign = fmpz_sgn(coeff);
    mp_limb_t cy;
    ulong limbs = (shift + bits) / FLINT_BITS;
    ulong rem_bits = (shift + bits) % FLINT_BITS;
    mp_limb_t mask;
    ulong size;

    if (sign == 0)  /* special case, deal with zero (store -borrow) */
    {
        if (borrow)
        {
            /* store -1 shifted and add save back in */
            arr[0] = ((~(mp_limb_t) 0) << shift) + save;

            /* com remaining limbs */
            if (limbs > 1)
                flint_mpn_store(arr + 1, limbs - 1, ~(mp_limb_t) 0);

            /* com remaining bits */
            if (limbs)
            {
                if (rem_bits)
                    arr[limbs] = (((mp_limb_t) 1) << rem_bits) - (mp_limb_t) 1;
            }
            else
            {
                /* mask off final limb */
                mask = (((mp_limb_t) 1) << rem_bits) - (mp_limb_t) 1;
                arr[limbs] &= mask;
            }

        }
        return borrow;
    }

    /*
       Let |c| = b. If c is -ve and negate == 0 or c is positive and negate is 1
       we want -b - borrow.

       If c is +ve and negate is 0 or c is negative and negate == 1, we want
       b - borrow.
     */
    if ((sign ^ negate) < 0)  /* -b - borrow = com(b) + 1 - borrow */
    {
        if (!COEFF_IS_MPZ(c))
        {
            /* compute d = -b - borrow */
            mp_limb_t d = (c < WORD(0) ? c - borrow : -c - borrow);

            /* store d << shift and add save back into place */
            arr[0] = (d << shift) + save;

            /* store carry from d<<shift and com remaining bits of second limb */
            if (limbs)
            {
                if (shift)
                    arr[1] =
                        (d >> (FLINT_BITS - shift)) +
                        ((~(mp_limb_t) 0) << shift);
                else
                    arr[1] = ~(mp_limb_t) 0;
            }

            size = 2;
        }
        else
        {
            __mpz_struct * ptr = COEFF_TO_PTR(c);
            size = FLINT_ABS(ptr->_mp_size);

            /* complement coefficient into arr */
            mpn_com(arr, ptr->_mp_d, size);

            /* deal with +1 - borrow */
            if (!borrow)
                mpn_add_1(arr, arr, size, 1);  /* cannot be a carry, else we com'd 0 */

            /* shift into place */
            if (shift)
            {
                cy = mpn_lshift(arr, arr, size, shift);
                if (limbs + (rem_bits != 0) > size)
                    arr[size++] = ((~(mp_limb_t) 0) << shift) + cy;
            }

            /* add back in saved bits from start of field */
            arr[0] += save;
        }

        if (limbs >= size)
        {
            /* com any additional limbs */
            if (limbs > size)
                flint_mpn_store(arr + size, limbs - size, ~(mp_limb_t) 0);

            /* com remaining bits */
            if (rem_bits)
                arr[limbs] = (((mp_limb_t) 1) << rem_bits) - (mp_limb_t) 1;
        }
        else
        {
            /* mask off final limb */
            mask = (((mp_limb_t) 1) << rem_bits) - (mp_limb_t) 1;
            arr[limbs] &= mask;
        }
        return 1;
    }
    else  /* b - borrow */
    {
        if (!COEFF_IS_MPZ(c))
        {
            /* compute d = b - borrow */
            mp_limb_t d = (c < WORD(0) ? -c - borrow : c - borrow);

            /* store d<<shift and add save back into place */
            arr[0] = (d << shift) + save;

            /* store carry from d<<shift */
            if (limbs + (rem_bits != 0) > 1)
            {
                if (shift)
                    arr[1] = (d >> (FLINT_BITS - shift));
            }
        }
        else
        {
            __mpz_struct *ptr = COEFF_TO_PTR(c);
            size = FLINT_ABS(ptr->_mp_size);

            /* shift into place */
            if (shift)
            {
                cy = mpn_lshift(arr, ptr->_mp_d, size, shift);
                if (cy)
                    arr[size++] = cy;
            }
            else
                flint_mpn_copyi(arr, ptr->_mp_d, size);

            /* deal with - borrow */
            if (borrow)
                mpn_sub_1(arr, arr, size, ((mp_limb_t) 1) << shift);

            /* add back in saved bits from start of field */
            arr[0] += save;
        }

        return 0;
    }
}
