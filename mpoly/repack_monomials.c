/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/*
    exps2 = packed vectors to unpack
    bits2 = bits per field in exps2
    len = number of packed vectors to unpack
    nfields = number of elements in each vector

    exps1 = destination for unpacked vector
    bits1 = number of bits desired in return packed vector

    return 1 if the repacking was successful, 0 if it failed
*/
int mpoly_repack_monomials(ulong * exps1, flint_bitcnt_t bits1,
                     const ulong * exps2, flint_bitcnt_t bits2, slong len,
                                                        const mpoly_ctx_t mctx)
{
    int success;
    slong i, j;
    slong nfields = mctx->nfields;
    slong N2 = mpoly_words_per_exp(bits2, mctx);
    slong N1 = mpoly_words_per_exp(bits1, mctx);
    TMP_INIT;

    if (bits1 == bits2)
    {
        for (i = 0; i < N2*len; i++)
            exps1[i] = exps2[i];
        return 1;
    }

    TMP_START;

    if (bits1 > bits2)
    {
        success = 1;

        if (bits1 <= FLINT_BITS && bits2 <= FLINT_BITS)
        {
            ulong * tmp_exps = (ulong *) TMP_ALLOC(nfields*sizeof(ulong));

            for (i = 0; i < len; i++)
            {
                mpoly_unpack_vec_ui(tmp_exps, exps2 + N2*i, bits2, nfields, 1);
                mpoly_pack_vec_ui(exps1 + N1*i, tmp_exps, bits1, nfields, 1);
            }
        }
        else
        {
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
    }
    else
    {
        success = 0;

        if (bits1 <= FLINT_BITS && bits2 <= FLINT_BITS)
        {
            ulong mask = (-UWORD(1)) << (bits1 - 1);
            ulong * tmp_exps = (ulong *) TMP_ALLOC(nfields*sizeof(ulong));
            for (i = 0; i < len; i++)
            {
                mpoly_unpack_vec_ui(tmp_exps, exps2 + N2*i, bits2, nfields, 1);
                for (j = 0; j < nfields; j++)
                {
                    if (tmp_exps[j] & mask)
                        goto cleanup;
                }
                mpoly_pack_vec_ui(exps1 + N1*i, tmp_exps, bits1, nfields, 1);
            }

            success = 1;
        }
        else
        {
            fmpz * tmp_exps = (fmpz *) TMP_ALLOC(nfields*sizeof(fmpz));

            for (i = 0; i < nfields; i++)
                fmpz_init(tmp_exps + i);

            for (i = 0; i < len; i++)
            {
                mpoly_unpack_vec_fmpz(tmp_exps, exps2 + N2*i, bits2, nfields, 1);
                if (_fmpz_vec_max_bits(tmp_exps, nfields) >= bits1)
                    goto cleanup1;
                mpoly_pack_vec_fmpz(exps1 + N1*i, tmp_exps, bits1, nfields, 1);
            }

            success = 1;
cleanup1:
            for (i = 0; i < nfields; i++)
                fmpz_clear(tmp_exps + i);
        }
    }

cleanup:

    TMP_END;

    return success;
}
