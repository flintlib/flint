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
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "perm.h"


mp_limb_t
_nmod_mat_det(nmod_mat_t A)
{
    mp_limb_t det;
    slong * P;

    slong m = A->r;
    slong rank;
    slong i;

    P = flint_malloc(sizeof(slong) * m);
    rank = nmod_mat_lu(P, A, 1);

    det = UWORD(0);

    if (rank == m)
    {
        det = UWORD(1);
        for (i = 0; i < m; i++)
            det = n_mulmod2_preinv(det, nmod_mat_entry(A, i, i),
                A->mod.n, A->mod.ninv);
    }

    if (_perm_parity(P, m) == 1)
        det = nmod_neg(det, A->mod);

    flint_free(P);
    return det;
}

mp_limb_t
nmod_mat_det(const nmod_mat_t A)
{
    nmod_mat_t tmp;
    mp_limb_t det;
    slong dim = A->r;

    if (dim != A->c)
    {
        flint_printf("Exception (nmod_mat_det). Non-square matrix.\n");
        flint_abort();
    }

    if (dim == 0) return UWORD(1);
    if (dim == 1) return nmod_mat_entry(A, 0, 0);

    nmod_mat_init_set(tmp, A);
    det = _nmod_mat_det(tmp);
    nmod_mat_clear(tmp);

    return det;
}
