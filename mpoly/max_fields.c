/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mpoly.h"


/* unpack the pointwise maximum of poly_exps into max_fields */
void mpoly_max_fields_ui(ulong * max_fields, const ulong * poly_exps,
                                          slong len, slong bits, slong nfields)
{
    slong i, N = words_per_exp(nfields, bits);
    ulong * pmax, mask;
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

    mpoly_unpack_vec(max_fields, pmax, bits, nfields, 1);

    TMP_END;
}


/*
    unpack the pointwise maximum of poly_exps into max_fields
    but write the results backwards (this backwards array is used for dense code)
*/
void mpoly_max_fields_ui_backwards(ulong * max_fields, const ulong * poly_exps,
                                          slong len, slong bits, slong nfields)
{
    slong i, j;

    mpoly_max_fields_ui(max_fields, poly_exps, len, bits, nfields);

    /* reverse the order */
    for (i = 0, j = nfields - 1; i < j; i++, j--)
    {
        ulong t = max_fields[j];
        max_fields[j] = max_fields[i];
        max_fields[i] = t;
    }
}

