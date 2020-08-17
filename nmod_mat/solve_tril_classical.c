/*
    Copyright (C) 2010,2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"
#include "nmod_vec.h"

void
nmod_mat_solve_tril_classical(nmod_mat_t X, const nmod_mat_t L,
                                                const nmod_mat_t B, int unit)
{
    int nlimbs;
    slong i, j, n, m;
    nmod_t mod;
    mp_ptr inv, tmp;

    n = L->r;
    m = B->c;
    mod = L->mod;

    if (!unit)
    {
        inv = _nmod_vec_init(n);
        for (i = 0; i < n; i++)
            inv[i] = n_invmod(nmod_mat_entry(L, i, i), mod.n);
    }
    else
        inv = NULL;

    nlimbs = _nmod_vec_dot_bound_limbs(n, mod);
    tmp = _nmod_vec_init(n);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
            tmp[j] = nmod_mat_entry(X, j, i);

        for (j = 0; j < n; j++)
        {
            mp_limb_t s;
            s = _nmod_vec_dot(L->rows[j], tmp, j, mod, nlimbs);
            s = nmod_sub(nmod_mat_entry(B, j, i), s, mod);
            if (!unit)
                s = n_mulmod2_preinv(s, inv[j], mod.n, mod.ninv);
            tmp[j] = s;
        }

        for (j = 0; j < n; j++)
            nmod_mat_entry(X, j, i) = tmp[j];
    }

    _nmod_vec_clear(tmp);
    if (!unit)
        _nmod_vec_clear(inv);
}
