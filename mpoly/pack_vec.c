/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* this file DOES NOT need to change with new orderings */

void mpoly_pack_vec_ui(ulong * exp1, const ulong * exp2, flint_bitcnt_t bits,
                                                      slong nfields, slong len)
{
    if (bits <= FLINT_BITS)
    {
        slong i, j;
        for (j = 0; j < len; j++)
        {
            ulong v = 0;
            slong shift = 0;
            i = 0;
            v |= *exp2++ << shift;
            shift += bits;      /* number of bits to encode 0th field */
            while (++i < nfields)
            {
                if (shift + bits > FLINT_BITS)
                {
                    *exp1++ = v;
                    v = 0;
                    shift = 0;
                }
                v |= *exp2++ << shift;
                shift += bits;      /* number of bits to encode ith field */
            }
            *exp1++ = v;
        }
    } else
    {
        slong j;
        ulong words_per_field = bits/FLINT_BITS;
        FLINT_ASSERT(bits%FLINT_BITS == 0);

        for (j = 0; j < len*nfields; j++, exp2++)
        {
            ulong size = 0;
            *exp1++ = *exp2;
            size++;
            while (size++ < words_per_field)
                *exp1++ = 0;
        }

    }
}

void mpoly_pack_vec_fmpz(ulong * exp1, const fmpz * exp2, flint_bitcnt_t bits,
                                                      slong nfields, slong len)
{
    if (bits <= FLINT_BITS)
    {
        slong i, j;
        for (j = 0; j < len; j++)
        {
            ulong v = 0;
            slong shift = 0;
            i = 0;
            FLINT_ASSERT(fmpz_abs_fits_ui(exp2));
            v |= fmpz_get_ui(exp2++) << shift;
            shift += bits;      /* number of bits to encode 0th field */
            while (++i < nfields)
            {
                if (shift + bits > FLINT_BITS)
                {
                    *exp1++ = v;
                    v = 0;
                    shift = 0;
                }
                FLINT_ASSERT(fmpz_abs_fits_ui(exp2));
                v |= fmpz_get_ui(exp2++) << shift;
                shift += bits;      /* number of bits to encode ith field */
            }
            *exp1++ = v;
        }

    } else
    {
        slong j;
        ulong words_per_field = bits/FLINT_BITS;
        FLINT_ASSERT(bits%FLINT_BITS == 0);

        for (j = 0; j < len*nfields; j++, exp2++)
        {
            ulong size = 0;
            if (fmpz_abs_fits_ui(exp2))
            {
                *exp1++ = fmpz_get_ui(exp2);
                size++;
            } else
            {
                __mpz_struct * mpz = COEFF_TO_PTR(*exp2);
                FLINT_ASSERT(mpz->_mp_size <= words_per_field);
                while (size < mpz->_mp_size)
                    *exp1++ = mpz->_mp_d[size++];
            }
            while (size++ < words_per_field)
                *exp1++ = 0;
        }
    }
}
