/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

void mpoly_pack_vec_ui(ulong * exp1, const ulong * exp2, slong bits, slong nfields, slong len) {
    slong i, j, shift;
    ulong v;
    for (j = 0; j < len; j++) {
        v = 0;
        shift = 0;
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
}

void mpoly_unpack_vec_ui(ulong * exp1, const ulong * exp2, slong bits, slong nfields, slong len) {
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
}
