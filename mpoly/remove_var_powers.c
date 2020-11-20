/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* this file does not need to change with new orderings */

/*
    set the array var_powers to the exponents on the monomial content of A
    and divide A by this monomial content
*/
void mpoly_remove_var_powers(
    fmpz * var_powers,
    ulong * Aexps,
    flint_bitcnt_t Abits,
    slong Alen,
    const mpoly_ctx_t mctx)
{
    slong i, N = mpoly_words_per_exp(Abits, mctx);
    fmpz * minfields;
    ulong * minexp;
    TMP_INIT;

    FLINT_ASSERT(Alen > 0);

    TMP_START;

    minexp = (ulong *) TMP_ALLOC(N*sizeof(ulong));
    minfields = (fmpz *) TMP_ALLOC(mctx->nfields*sizeof(fmpz));

    for (i = 0; i < mctx->nfields; i++)
        fmpz_init(minfields + i);
    mpoly_min_fields_fmpz(minfields, Aexps, Alen, Abits, mctx);

    mpoly_get_monomial_ffmpz_unpacked_ffmpz(var_powers, minfields, mctx);

    mpoly_set_monomial_ffmpz(minexp, var_powers, Abits, mctx);

    if (!mpoly_monomial_is_zero(minexp, N))
    {
        if (Abits <= FLINT_BITS)
        {
            for (i = 0; i < Alen; i++)
                mpoly_monomial_sub(Aexps + N*i, Aexps + N*i, minexp, N);
        }
        else
        {
            for (i = 0; i < Alen; i++)
                mpoly_monomial_sub_mp(Aexps + N*i, Aexps + N*i, minexp, N);
        }
    }

    TMP_END;
}

