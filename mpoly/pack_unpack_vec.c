/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

void mpoly_pack_vec(ulong * exp1, const ulong * exp2, slong bits, slong nfields, slong len) {
    slong i, j, shift, fields_per_word = FLINT_BITS/bits;
    ulong v;
    for (j = 0; j < len; j++) {
        v = 0;
        shift = bits*fields_per_word;
        for (i = 0; i < nfields; i++) {
            shift -= bits;
            v |= *exp2++ << shift;
            if (shift == 0)
            {
                *exp1++ = v;
                v = 0;
                shift = bits*fields_per_word;
            }
        }
        if (shift != bits*fields_per_word)
        {
            *exp1++ = v;
        }
    }
}

void mpoly_unpack_vec(ulong * exp1, const ulong * exp2, slong bits, slong nfields, slong len) {
    slong i, j, shift, fields_per_word = FLINT_BITS/bits;
    ulong u, mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    for (j = 0; j < len; j++) {
        u = *exp2++;
        shift = bits*fields_per_word;
        for (i = 0; i < nfields; i++) {
            shift -= bits;
            *exp1++ = (u >> shift) & mask;
            if (shift == 0 && i + 1 < nfields)
            {
                u = *exp2++;
                shift = bits*fields_per_word;
            }
        }
    }
}
