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

#include "flint-mparam.h"
#include "ulong_extras.h"
#include "nmod_mat.h"
#include "fmpz.h"

void
nmod_mat_scalar_addmul_ui(nmod_mat_t C, const nmod_mat_t A,
                                          const nmod_mat_t B, const ulong c)
{
    slong i, j;

    if (c == UWORD(0))
    {
        if (C != A)
            nmod_mat_set(C, A);
        return;
    }

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < A->c; j++)
        {
            nmod_mat_entry(C, i, j) =
                n_addmod(nmod_mat_entry(A, i, j),
                         n_mulmod2_preinv(nmod_mat_entry(B, i, j), c,
                                          B->mod.n, B->mod.ninv),
                         A->mod.n);
        }
    }
}

void
nmod_mat_scalar_mul(nmod_mat_t B, const nmod_mat_t A, ulong c)
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
    else if (A->r * A->c > FLINT_MULMOD_SHOUP_THRESHOLD && A->mod.norm > 0)
    {
        slong i, j;
        ulong c_pr = n_mulmod_precomp_shoup(c, A->mod.n);

        for (i = 0; i < A->r; i++)
            for (j = 0; j < A->c; j++)
                nmod_mat_entry(B, i, j) = n_mulmod_shoup(
                    c, nmod_mat_entry(A, i, j), c_pr, A->mod.n);
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
nmod_mat_scalar_mul_fmpz(nmod_mat_t B, const nmod_mat_t A, const fmpz_t c)
{
    nmod_mat_scalar_mul(B, A, fmpz_get_nmod(c, B->mod));
}
