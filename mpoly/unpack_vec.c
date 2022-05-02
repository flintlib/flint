/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

void mpoly_unpack_vec_ui(ulong * exp1, const ulong * exp2, flint_bitcnt_t bits,
                                                     slong nfields, slong len)
{
    if (bits < FLINT_BITS)
    {
        slong i, j, shift;
        ulong u, mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        for (j = 0; j < len; j++)
        {
            i = 0;
            u = *exp2++;
            shift = 0;
            *exp1++ = u & mask;
            u = u >> bits;      /* number of bits to encode 0th field */
            shift += bits;      /* number of bits to encode 0th field */
            while (++i < nfields)
            {
                if (shift + bits > FLINT_BITS) 
                {
                    u = *exp2++;
                    shift = 0;
                }
                *exp1++ = u & mask;
                u = u >> bits;      /* number of bits to encode ith field */
                shift += bits;      /* number of bits to encode ith field */
            }
        }
    }
    else
    {
        slong j;
        ulong words_per_field = bits/FLINT_BITS;
        FLINT_ASSERT(bits%FLINT_BITS == 0);

        for (j = 0; j < len*nfields; j++, exp2 += words_per_field)
        {
            *exp1++ = *exp2;
        }
    }
}

void mpoly_unpack_vec_fmpz(fmpz * exp1, const ulong * exp2, flint_bitcnt_t bits,
                                                      slong nfields, slong len)
{
    if (bits < FLINT_BITS)
    {
        slong i, j, shift;
        ulong u, mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        for (j = 0; j < len; j++)
        {
            i = 0;
            u = *exp2++;
            shift = 0;
            fmpz_set_ui(exp1++, u & mask);
            u = u >> bits;      /* number of bits to encode 0th field */
            shift += bits;      /* number of bits to encode 0th field */
            while (++i < nfields)
            {
                if (shift + bits > FLINT_BITS)
                {
                    u = *exp2++;
                    shift = 0;
                }
                fmpz_set_ui(exp1++, u & mask);
                u = u >> bits;      /* number of bits to encode ith field */
                shift += bits;      /* number of bits to encode ith field */
            }
        }
    }
    else
    {
        slong j;
        ulong words_per_field = bits/FLINT_BITS;
        FLINT_ASSERT(bits%FLINT_BITS == 0);

        for (j = 0; j < len*nfields; j++, exp2 += words_per_field)
        {
            ulong size = words_per_field;
            while (size > 1 && exp2[size - 1] == 0)
                size--;
            if (size == 1)
            {
                fmpz_set_ui(exp1, exp2[0]);
            }
            else
            {
                __mpz_struct * mpz = _fmpz_promote(exp1);
                if (mpz->_mp_alloc < words_per_field)
                    mpz_realloc2(mpz, bits);
                mpz->_mp_size = size;
                flint_mpn_copyi(mpz->_mp_d, exp2, size);
            }
            exp1++;
        }
    }
}
