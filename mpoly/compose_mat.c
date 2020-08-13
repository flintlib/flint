/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* unfortunate function missing from fmpz_mat */
void fmpz_mat_mul_vec(fmpz * v, const fmpz_mat_t M, fmpz * u)
{
    slong i;
    slong r = fmpz_mat_nrows(M);
    slong c = fmpz_mat_ncols(M);

    for (i = 0; i < r; i++)
        _fmpz_vec_dot(v + i, M->rows[i], u, c);
}

/* Fill the compose matrix for the T_mpoly_compose_T_mpoly_gen functions */
void mpoly_compose_mat_gen(fmpz_mat_t M, const slong * c,
                             const mpoly_ctx_t mctxB, const mpoly_ctx_t mctxAC)
{
    slong i, gi, j;
    fmpz * t;

    FLINT_ASSERT(fmpz_mat_nrows(M) == mctxAC->nfields + 1);
    FLINT_ASSERT(fmpz_mat_ncols(M) == mctxB->nfields);

    fmpz_mat_zero(M);

    t = _fmpz_vec_init(mctxAC->nfields);

    for (i = 0; i < mctxB->nvars; i++)
    {
        gi = mpoly_gen_index(i, mctxB);
        if (0 <= c[i] && c[i] < mctxAC->nfields)
        {
            mpoly_gen_fields_fmpz(t, c[i], mctxAC);
            for (j = 0; j < mctxAC->nfields; j++)
                fmpz_swap(fmpz_mat_entry(M, j, gi), t + j);
        }
        else
        {
            /* c[i] corresponds to zero */
            fmpz_one(fmpz_mat_entry(M, mctxAC->nfields, gi));
        }
    }

    _fmpz_vec_clear(t, mctxAC->nfields);
}

/*
    Fill the column of the compose matrix for the variable Bvar assuming
    Bvar maps to the monomial in (Cexps, Cbits).
    NULL may be passed in for Cexp to indicate that Bbar maps to 0.
*/
void mpoly_compose_mat_fill_column(fmpz_mat_t M, const ulong * Cexp,
                     flint_bitcnt_t Cbits, slong Bvar, const mpoly_ctx_t mctxB,
                                                      const mpoly_ctx_t mctxAC)
{
    slong Bidx = mpoly_gen_index(Bvar, mctxB);
    slong j;
    fmpz * t;

    if (Cexp == NULL)
    {
        j = mctxAC->nfields;
        fmpz_one(fmpz_mat_entry(M, j, Bidx));

        for (j--; j >= 0; j--)
            fmpz_zero(fmpz_mat_entry(M, j, Bidx));

        return;
    }

    t = _fmpz_vec_init(mctxAC->nfields);

    mpoly_unpack_vec_fmpz(t, Cexp, Cbits, mctxAC->nfields, 1);

    j = mctxAC->nfields;
    fmpz_zero(fmpz_mat_entry(M, j, Bidx));
    for (j--; j >= 0; j--)
        fmpz_swap(fmpz_mat_entry(M, j, Bidx), t + j);

    _fmpz_vec_clear(t, mctxAC->nfields);
}
