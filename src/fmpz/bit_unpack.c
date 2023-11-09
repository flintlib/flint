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
fmpz_bit_unpack(fmpz_t coeff, mp_srcptr arr, flint_bitcnt_t shift,
                flint_bitcnt_t bits, int negate, int borrow)
{
    mp_limb_t mask, sign;
    ulong limbs = (shift + bits) / FLINT_BITS;
    ulong rem_bits = (shift + bits) % FLINT_BITS;

    /* determine if field is positive or negative */
    if (rem_bits)
        sign = ((((mp_limb_t) 1) << (rem_bits - 1)) & arr[limbs]);
    else
        sign = ((((mp_limb_t) 1) << (FLINT_BITS - 1)) & arr[limbs - 1]);

    if (bits <= SMALL_FMPZ_BITCOUNT_MAX)  /* fits into a small coeff */
    {
        _fmpz_demote(coeff);

        /* mask for the given number of bits */
        mask = (((mp_limb_t) 1) << bits) - (mp_limb_t) 1;

        if (limbs + (rem_bits != 0) > 1)  /* field crosses a limb boundary */
            (*coeff) =
                ((arr[0] >> shift) + (arr[1] << (FLINT_BITS - shift))) & mask;
        else  /* field is in first limb only, mask it */
            (*coeff) = (arr[0] >> shift) & mask;

        /* sign extend */
        if (sign)
            (*coeff) += ((~(mp_limb_t) 0) << bits);

        /* determine whether we need to return a borrow */
        sign = (*coeff < (mp_limb_signed_t) 0 ? (mp_limb_t) 1 : (mp_limb_t) 0);

        /* deal with borrow */
        if (borrow)
        {
            (*coeff)++;
            if ((*coeff) > COEFF_MAX)
            {
                fmpz v = *coeff;
                *coeff = 0;
                fmpz_set_ui(coeff, v);
            }
        }

        /* negate the coeff if necessary */
        if (negate)
            fmpz_neg(coeff, coeff);

        return (sign != (mp_limb_t) 0);
    }
    else  /* large coefficient */
    {
        __mpz_struct * mcoeff;
        mp_limb_t * p;
        ulong l, b;

        mcoeff = _fmpz_promote(coeff);

        /* the number of limbs to hold the bitfield _including_ b extra bits */
        l = (bits - 1) / FLINT_BITS + 1;
        b = bits % FLINT_BITS;

        /* allocate space for l limbs only */
        mpz_realloc(mcoeff, l);
        p = mcoeff->_mp_d;

        /* shift in l limbs */
        if (shift)
            mpn_rshift(p, arr, l, shift);
        else
            flint_mpn_copyi(p, arr, l);

        /* shift in any rem_bits that weren't already shifted */
        if (limbs + (rem_bits != 0) > l)
            p[l - 1] += (arr[limbs] << (FLINT_BITS - shift));

        /* mask off the last limb, if not full */
        if (b)
        {
            mask = (((mp_limb_t) 1) << b) - (mp_limb_t) 1;
            p[l - 1] &= mask;
        }

        if (sign != (mp_limb_t) 0)
        {
            /* sign extend */
            if (b)
                p[l - 1] += ((~(mp_limb_t) 0) << b);

            /* negate */
            mpn_com(p, p, l);
            if (!borrow)
                mpn_add_1(p, p, l, 1);

            /* normalise */
            while (l && (p[l - 1] == (mp_limb_t) 0))
                l--;

            mcoeff->_mp_size = -l;

            sign = 1;
        }
        else
        {
            /* deal with borrow */
            if (borrow)
                mpn_add_1(p, p, l, 1);

            /* normalise */
            while (l && (p[l - 1] == (mp_limb_t) 0))
                l--;

            mcoeff->_mp_size = l;
            sign = 0;
        }

        /* negate if required */
        if (negate)
            mpz_neg(mcoeff, mcoeff);

        /* coeff may fit in a small */
        _fmpz_demote_val(coeff);

        return sign;
    }
}

void
fmpz_bit_unpack_unsigned(fmpz_t coeff, mp_srcptr arr,
                         flint_bitcnt_t shift, flint_bitcnt_t bits)
{
    ulong limbs = (shift + bits) / FLINT_BITS;
    ulong rem_bits = (shift + bits) % FLINT_BITS;
    mp_limb_t mask;

    if (bits <= SMALL_FMPZ_BITCOUNT_MAX)  /* fits into a small coeff */
    {
        _fmpz_demote(coeff);

        /* mask for the given number of bits */
        mask = (((mp_limb_t) 1) << bits) - (mp_limb_t) 1;

        if (limbs + (rem_bits != 0) > 1)  /* field crosses a limb boundary */
            (*coeff) =
                ((arr[0] >> shift) + (arr[1] << (FLINT_BITS - shift))) & mask;
        else  /* field is in first limb only, mask it */
            (*coeff) = (arr[0] >> shift) & mask;
    }
    else  /* large coefficient */
    {
        __mpz_struct * mcoeff;
        mp_limb_t * p;
        ulong l, b;

        mcoeff = _fmpz_promote(coeff);

        /* the number of limbs to hold the bitfield _including_ b extra bits */
        l = (bits - 1) / FLINT_BITS + 1;
        b = bits % FLINT_BITS;

        /* allocate space for l limbs only */
        mpz_realloc(mcoeff, l);
        p = mcoeff->_mp_d;

        /* shift in l limbs */
        if (shift)
            mpn_rshift(p, arr, l, shift);
        else
            flint_mpn_copyi(p, arr, l);

        /* shift in any rem_bits that weren't already shifted */
        if (limbs + (rem_bits != 0) > l)
            p[l - 1] += (arr[limbs] << (FLINT_BITS - shift));

        /* mask off the last limb, if not full */
        if (b)
        {
            mask = (((mp_limb_t) 1) << b) - (mp_limb_t) 1;
            p[l - 1] &= mask;
        }

        /* normalise */
        while (l && (p[l - 1] == (mp_limb_t) 0))
            l--;

        mcoeff->_mp_size = l;

        /* coeff may fit in a small */
        _fmpz_demote_val(coeff);
    }
}
