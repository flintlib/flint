#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "mpoly.h"

/* maximize "bits" while keeping "(nfields-1)/(FLINT_BITS/bits)+1" constant */
slong mpoly_optimize_bits(slong bits, slong nfields) {
    while (bits < FLINT_BITS &&   (nfields - 1)/(FLINT_BITS/(bits    ))
                               == (nfields - 1)/(FLINT_BITS/(bits + 1))
          )
    {
        bits++;
    }

    return bits;
}



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

/*
    compute number of bits required to store user_exp in packed format
    the returned number of bits includes space for a zero'd signed bit
    a return value of > FLINT_BITS indicates an error (it doesn't fit)
*/
slong mpoly_exp_bits(const ulong * user_exp, slong nfields, int deg)
{
    slong i, bits, exp_bits = 8;
    ulong max = 0;
    if (deg)
    {
        for (i = 0; i < nfields - 1; i++)
        {
            max += user_exp[i];
            if (max < user_exp[i])
                return FLINT_BITS + 1;
        }
    } else
    {
        for (i = 0; i < nfields; i++)
        {
            if (max < user_exp[i])
                max = user_exp[i];
        }
    }

    bits = FLINT_BIT_COUNT(max);
    while (bits >= exp_bits)
        exp_bits += 1;

    return exp_bits;
}
