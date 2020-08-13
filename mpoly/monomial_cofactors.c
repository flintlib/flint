/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/*
    If monomials can be factored of A and B leaving identical exponents,
    store the monomials in Abarexps and Bbarsexps and return 1. Otherwise,
    return 0.
*/
int mpoly_monomial_cofactors(
    fmpz * Abarexps,
    fmpz * Bbarexps,
    const ulong * Aexps, flint_bitcnt_t Abits,
    const ulong * Bexps, flint_bitcnt_t Bbits,
    slong length,
    const mpoly_ctx_t mctx)
{
    int success;
    slong i, j;
    slong nvars = mctx->nvars;
    slong NA = mpoly_words_per_exp(Abits, mctx);
    slong NB = mpoly_words_per_exp(Bbits, mctx);
    fmpz * Aexp, * Bexp, * minAexp, * minBexp;
    fmpz_t t1, t2;
    TMP_INIT;

    FLINT_ASSERT(length > 0);

    fmpz_init(t1);
    fmpz_init(t2);

    TMP_START;

    Aexp = (fmpz *) TMP_ALLOC(4*nvars*sizeof(fmpz));
    Bexp = Aexp + 1*nvars;
    minAexp = Aexp + 2*nvars;
    minBexp = Aexp + 3*nvars;

    for (j = 0; j < nvars; j++)
    {
        fmpz_init(Aexp + j);
        fmpz_init(Bexp + j);
        fmpz_init(minAexp + j);
        fmpz_init(minBexp + j);
    }

    /* {A|B}barexps holds the leading exponents */
    mpoly_get_monomial_ffmpz(Abarexps, Aexps + NA*0, Abits, mctx);
    mpoly_get_monomial_ffmpz(Bbarexps, Bexps + NB*0, Bbits, mctx);
    _fmpz_vec_set(minAexp, Abarexps, nvars);
    _fmpz_vec_set(minBexp, Bbarexps, nvars);

    for (i = 0; i < length; i++)
    {
        mpoly_get_monomial_ffmpz(Aexp, Aexps + NA*i, Abits, mctx);
        mpoly_get_monomial_ffmpz(Bexp, Bexps + NB*i, Bbits, mctx);
        _fmpz_vec_min_inplace(minAexp, Aexp, nvars);
        _fmpz_vec_min_inplace(minBexp, Bexp, nvars);
        for (j = 0; j < nvars; j++)
        {
            fmpz_add(t1, Abarexps + j, Bexp + j);
            fmpz_add(t2, Bbarexps + j, Aexp + j);
            success = fmpz_equal(t1, t2);
            if (!success)
                goto cleanup;
        }
    }

    /* put {A|B}'s cofactor monomial in {A|B}barexps */
    _fmpz_vec_max(Bbarexps, minAexp, minBexp, nvars);
    _fmpz_vec_sub(Abarexps, Bbarexps, minBexp, nvars);
    _fmpz_vec_sub(Bbarexps, Bbarexps, minAexp, nvars);

    success = 1;

cleanup:

    for (j = 0; j < nvars; j++)
    {
        fmpz_clear(Aexp + j);
        fmpz_clear(Bexp + j);
        fmpz_clear(minAexp + j);
        fmpz_clear(minBexp + j);
    }

    TMP_END;

    fmpz_clear(t1);
    fmpz_clear(t2);

    return success;
}
