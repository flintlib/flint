/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/*
    exps2 = packed vectors to unpack
    bits2 = bits per field in exps2
    len = number of packed vectors to unpack
    nfields = number of elements in each vector

    exps1 = destination for unpacked vector
    bits1 = number of bits desired in return packed vector
*/
void mpoly_repack_monomials(ulong * exps1, slong bits1,
                      const ulong * exps2, slong bits2, slong len,
                                                        const mpoly_ctx_t mctx)
{
    slong i;
    slong nfields = mctx->nfields;
    slong N2 = mpoly_words_per_exp(bits2, mctx);
    slong N1 = mpoly_words_per_exp(bits1, mctx);
    TMP_INIT;

    if (bits1 == bits2) {
        for (i = 0; i < N2*len; i++)
            exps1[i] = exps2[i];
        return;
    }

    TMP_START;

    if (bits1 <= FLINT_BITS && bits2 <= FLINT_BITS)
    {
        ulong * tmp_exps = (ulong *) TMP_ALLOC(nfields*sizeof(ulong));

        for (i = 0; i < len; i++)
        {
            mpoly_unpack_vec_ui(tmp_exps, exps2 + N2*i, bits2, nfields, 1);
            mpoly_pack_vec_ui(exps1 + N1*i, tmp_exps, bits1, nfields, 1);
        }

    } else {

        fmpz * tmp_exps = (fmpz *) TMP_ALLOC(nfields*sizeof(fmpz));

        for (i = 0; i < nfields; i++)
            fmpz_init(tmp_exps + i);

        for (i = 0; i < len; i++)
        {
            mpoly_unpack_vec_fmpz(tmp_exps, exps2 + N2*i, bits2, nfields, 1);
            mpoly_pack_vec_fmpz(exps1 + N1*i, tmp_exps, bits1, nfields, 1);
        }

        for (i = 0; i < nfields; i++)
            fmpz_clear(tmp_exps + i);

    }

    TMP_END;
}
