/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "assert.h"

void mpoly_pack_vec_ui(ulong * exp1, const ulong * exp2, slong bits, slong nfields, slong len) {

    if (bits <= FLINT_BITS)
    {
        slong i, j;
        for (j = 0; j < len; j++) {
            ulong v = 0;
            slong shift = 0;
            i = 0;
            v |= *exp2++ << shift;
            shift += bits;      /* number of bits to encode 0th field */
            while (++i < nfields) {
                if (shift + bits > FLINT_BITS) {
                    *exp1++ = v;
                    v = 0;
                    shift = 0;
                }
                v |= *exp2++ << shift;
                shift += bits;      /* number of bits to encode ith field */
            }
            *exp1++ = v;
        }

    } else {

        slong j;
        ulong words_per_field = bits/FLINT_BITS;
        assert(bits%FLINT_BITS == 0);

        for (j = 0; j < len*nfields; j++, exp2++) {

                ulong size = 0;
                *exp1++ = *exp2;
                size++;
                while (size++ < words_per_field)
                    *exp1++ = 0;
        }

    }
}

void mpoly_pack_vec_fmpz(ulong * exp1, const fmpz * exp2, mp_bitcnt_t bits, slong nfields, slong len) {

    if (bits <= FLINT_BITS)
    {
        slong i, j;
        for (j = 0; j < len; j++) {
            ulong v = 0;
            slong shift = 0;
            i = 0;
            assert(fmpz_abs_fits_ui(exp2));
            v |= fmpz_get_ui(exp2++) << shift;
            shift += bits;      /* number of bits to encode 0th field */
            while (++i < nfields) {
                if (shift + bits > FLINT_BITS) {
                    *exp1++ = v;
                    v = 0;
                    shift = 0;
                }
                assert(fmpz_abs_fits_ui(exp2));
                v |= fmpz_get_ui(exp2++) << shift;
                shift += bits;      /* number of bits to encode ith field */
            }
            *exp1++ = v;
        }

    } else {

        slong j;
        ulong words_per_field = bits/FLINT_BITS;
        assert(bits%FLINT_BITS == 0);

        for (j = 0; j < len*nfields; j++, exp2++) {

                ulong size = 0;
                if (fmpz_abs_fits_ui(exp2)) {
                    *exp1++ = fmpz_get_ui(exp2);
                    size++;
                } else {
                    __mpz_struct * mpz = COEFF_TO_PTR(*exp2);
                    assert(mpz->_mp_size <= words_per_field);
                    while (size < mpz->_mp_size)
                        *exp1++ = mpz->_mp_d[size++];
                }
                while (size++ < words_per_field)
                    *exp1++ = 0;
        }
    }
}

void mpoly_unpack_vec_ui(ulong * exp1, const ulong * exp2, slong bits, slong nfields, slong len) {

    if (bits <= FLINT_BITS) {

        slong i, j, shift;
        ulong u, mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        for (j = 0; j < len; j++) {
            i = 0;
            u = *exp2++;
            shift = 0;
            *exp1++ = u & mask;
            u = u >> bits;      /* number of bits to encode 0th field */
            shift += bits;      /* number of bits to encode 0th field */
            while (++i < nfields) {
                if (shift + bits > FLINT_BITS) {
                    u = *exp2++;
                    shift = 0;
                }
                *exp1++ = u & mask;
                u = u >> bits;      /* number of bits to encode ith field */
                shift += bits;      /* number of bits to encode ith field */
            }
        }

    } else {

        slong j;
        ulong words_per_field = bits/FLINT_BITS;
        assert(bits%FLINT_BITS == 0);

        for (j = 0; j < len*nfields; j++, exp2 += words_per_field) {
                *exp1++ = *exp2;
        }

    }
}

void mpoly_unpack_vec_fmpz(fmpz * exp1, const ulong * exp2, mp_bitcnt_t bits, slong nfields, slong len) {

    if (bits <= FLINT_BITS) {

        slong i, j, shift;
        ulong u, mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        for (j = 0; j < len; j++) {
            i = 0;
            u = *exp2++;
            shift = 0;
            fmpz_set_ui(exp1++, u & mask);
            u = u >> bits;      /* number of bits to encode 0th field */
            shift += bits;      /* number of bits to encode 0th field */
            while (++i < nfields) {
                if (shift + bits > FLINT_BITS) {
                    u = *exp2++;
                    shift = 0;
                }
                fmpz_set_ui(exp1++, u & mask);
                u = u >> bits;      /* number of bits to encode ith field */
                shift += bits;      /* number of bits to encode ith field */
            }
        }

    } else {

        slong j;
        ulong words_per_field = bits/FLINT_BITS;
        assert(bits%FLINT_BITS == 0);


        for (j = 0; j < len*nfields; j++, exp2 += words_per_field) {

                ulong size = words_per_field;
                while (size > 1 && exp2[size - 1] == 0)
                    size--;
                if (size == 1) {
                    fmpz_set_ui(exp1, exp2[0]);
                } else {
                    __mpz_struct * mpz = _fmpz_promote(exp1);
                    mpz_realloc2(mpz, bits);
                    mpn_copyi(mpz->_mp_d, exp2, size);
                    mpz->_mp_size = size;
                }

                exp1++;
        }
    }
}
