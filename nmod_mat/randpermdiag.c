/*
    Copyright (C) 2010,2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "ulong_extras.h"
#include "nmod_mat.h"
#include "nmod_vec.h"
#include "perm.h"

int
nmod_mat_randpermdiag(nmod_mat_t mat, flint_rand_t state,
                            mp_srcptr diag, slong n)
{
    int parity;
    slong i;
    slong * rows;
    slong * cols;

    rows = _perm_init(mat->r);
    cols = _perm_init(mat->c);

    parity = _perm_randtest(rows, mat->r, state);
    parity ^= _perm_randtest(cols, mat->c, state);

    nmod_mat_zero(mat);
    for (i = 0; i < n; i++)
        nmod_mat_entry(mat, rows[i], cols[i]) = diag[i];

    _perm_clear(rows);
    _perm_clear(cols);

    return parity;
}
