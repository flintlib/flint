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
    slong i, j, N = words_per_exp(nfields, bits);
    ulong * pmax, mask, t;
    TMP_INIT;

    mask = 0;
    for (i = 0; i < FLINT_BITS/bits; i++)
        mask = (mask << bits) + (UWORD(1) << (bits - 1));

    TMP_START;

    pmax = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    for (i = 0; i < N; i++)
        pmax[i] = 0;
    for (i = 0; i < len; i++)
        mpoly_monomial_max(pmax, pmax, poly_exps + i*N, bits, N, mask);

    mpoly_unpack_vec(max_degs, pmax, bits, nfields, 1);

    /* reverse the order */
    for (i = 0, j = nfields - 1; i < j; i++, j--)
    {
        t = max_degs[j];
        max_degs[j] = max_degs[i];
        max_degs[i] = t;
    }

    TMP_END;
}
