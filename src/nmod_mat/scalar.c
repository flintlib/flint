/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Ashish Kedia
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_mat.h"
#include "fmpz.h"

void
nmod_mat_scalar_addmul_ui(nmod_mat_t dest, const nmod_mat_t X,
                                          const nmod_mat_t Y, const mp_limb_t b)
{
    slong i, j;

    if (b == UWORD(0))
    {
        if (dest != X)
            nmod_mat_set(dest, X);
        return;
    }

    for (i = 0; i < X->r; i++)
    {
        for (j = 0; j < X->c; j++)
        {
            nmod_mat_entry(dest, i, j) =
                n_addmod(nmod_mat_entry(X, i, j),
                    n_mulmod2_preinv(nmod_mat_entry(Y, i, j), b, Y->mod.n,
                                                        Y->mod.ninv), X->mod.n);
        }
    }
}

#define UWORD_HALF (UWORD_MAX / 2 + 1)

void
nmod_mat_scalar_mul(nmod_mat_t B, const nmod_mat_t A, mp_limb_t c)
{
    if (c == UWORD(0))
    {
        nmod_mat_zero(B);
    }
    else if (c == UWORD(1))
    {
        nmod_mat_set(B, A);
    }
    else if (c == A->mod.n - UWORD(1))
    {
        nmod_mat_neg(B, A);
    }
    else if (A->r * A->c > 10 && A->mod.n < UWORD_HALF)
    {
        slong i, j;
        mp_limb_t w_pr = n_mulmod_precomp_shoup(c, A->mod.n);

        for (i = 0; i < A->r; i++)
            for (j = 0; j < A->c; j++)
                nmod_mat_entry(B, i, j) = n_mulmod_shoup(
                    c, nmod_mat_entry(A, i, j), w_pr, A->mod.n);
    }
    else
    {
        slong i, j;

        for (i = 0; i < A->r; i++)
            for (j = 0; j < A->c; j++)
                nmod_mat_entry(B, i, j) = n_mulmod2_preinv(
                    nmod_mat_entry(A, i, j), c, A->mod.n, A->mod.ninv);
    }
}

void
nmod_mat_scalar_mul_fmpz(nmod_mat_t res, const nmod_mat_t M, const fmpz_t c)
{
    nmod_mat_scalar_mul(res, M, fmpz_get_nmod(c, res->mod));
}
