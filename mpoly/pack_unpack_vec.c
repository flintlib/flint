#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "mpoly.h"
#include <assert.h>


void mpoly_pack_vec(ulong * exp1, const ulong * exp2, slong bits, slong fields, slong len) {
    slong i, j, shift, fields_per_word = FLINT_BITS/bits;
    ulong v, mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    for (j = 0; j < len; j++) {
        v = 0;
        shift = bits*fields_per_word;
        for (i = 0; i < fields; i++) {
            shift -= bits;
            assert((*exp2 & mask) == *exp2);
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

void mpoly_unpack_vec(ulong * exp1, const ulong * exp2, slong bits, slong fields, slong len) {
    slong i, j, shift, fields_per_word = FLINT_BITS/bits;
    ulong u, mask = (-UWORD(1)) >> (FLINT_BITS - bits);
    for (j = 0; j < len; j++) {
        u = *exp2++;
        shift = bits*fields_per_word;
        for (i = 0; i < fields; i++) {
            shift -= bits;
            *exp1++ = (u >> shift) & mask;
            if (shift == 0 && i + 1 < fields)
            {
                u = *exp2++;
                shift = bits*fields_per_word;
            }
        }
    }
}
