/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "mpoly.h"

/*
    exps2 = packed vectors to unpack
    bits2 = bits per field in exps2
    len = number of packed vectors to unpack
    nfields = number of elements in each vector

    exps1 = destination for unpacked vector
    bits1 = number of bits desired in return packed vector
*/
void mpoly_unpack_monomials(ulong * exps1, slong bits1,
                   const ulong * exps2, slong bits2, slong len, slong nfields)
{
    slong i, N2 = words_per_exp(nfields, bits2);
    slong    N1 = words_per_exp(nfields, bits1);
    ulong * tmp_exps;
    TMP_INIT;

    FLINT_ASSERT(bits1 >= bits2);
    if (bits1 == bits2) {
        for (i = 0; i < N2*len; i++)
            exps1[i] = exps2[i];
        return;
    }

    TMP_START;
    tmp_exps = (ulong *) TMP_ALLOC(nfields*sizeof(ulong));

    for (i = 0; i < len; i++)
    {
        mpoly_unpack_vec(tmp_exps, exps2 + i*N2, bits2, nfields, 1);
        mpoly_pack_vec(exps1 + i*N1, tmp_exps, bits1, nfields, 1);
    }

    TMP_END;
}
