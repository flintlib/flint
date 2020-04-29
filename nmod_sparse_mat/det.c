/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_sparse_mat.h"
#include "perm.h"

mp_limb_t
nmod_sparse_mat_det(const nmod_sparse_mat_t M)
{
    slong rank, i;
    slong *P, *Q;
    mp_limb_t det;
    nmod_sparse_mat_t L, U;
    if (M->r != M->c)
    {
        flint_printf("Exception (nmod_mat_det). Non-square matrix.\n");
        flint_abort();
    }

    if (M->r == 0) return UWORD(1);
    if (nmod_sparse_mat_is_zero(M)) return UWORD(0);
    if (M->r == 1) return M->rows[0].entries[0].val;


    P = flint_malloc(M->r*sizeof(*P));
    Q = flint_malloc(M->c*sizeof(*P));
    nmod_sparse_mat_init(L, M->r, M->c, M->mod);
    nmod_sparse_mat_init(U, M->r, M->c, M->mod);
    rank = nmod_sparse_mat_lu(P, Q, L, U, M);

    det = UWORD(0);

    if (rank == M->r)
    {
        det = UWORD(1);
        for (i = 0; i < M->r; i++)
            det = n_mulmod2_preinv(det, U->rows[i].entries[0].val, M->mod.n, M->mod.ninv);
        if ((_perm_parity(P, M->r) == 1) ^ (_perm_parity(Q, M->c) == 1))
            det = nmod_neg(det, M->mod);
    }

    flint_free(P);
    flint_free(Q);
    nmod_sparse_mat_clear(L);
    nmod_sparse_mat_clear(U);

    return det;
}
