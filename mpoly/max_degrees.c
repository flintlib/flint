/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"


void mpoly_max_degrees(ulong * max_degs, const ulong * poly_exps,
                                          slong len, slong bits, slong nfields)
{
    slong i, j, k, N = words_per_exp(nfields, bits);
    ulong * tmp_exps;
    TMP_INIT;

    for (i = 0; i < nfields; i++)
        max_degs[i] = 0;

    TMP_START;
    tmp_exps = (ulong *) TMP_ALLOC(nfields*sizeof(ulong));

    for (i = 0; i < len; i++)
    {
        mpoly_unpack_vec(tmp_exps, poly_exps + i*N, bits, nfields, 1);
        for (j = 0, k = nfields - 1; j < nfields; j++, k--)
        {
            if (max_degs[j] < tmp_exps[k])
                max_degs[j] = tmp_exps[k];
        }
    }

    TMP_END;
}
