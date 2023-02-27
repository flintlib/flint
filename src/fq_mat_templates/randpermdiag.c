/*
    Copyright (C) 2010,2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include "flint.h"
#include "ulong_extras.h"
#include "perm.h"

int
TEMPLATE(T, mat_randpermdiag) (TEMPLATE(T, mat_t) mat, flint_rand_t state,
                               TEMPLATE(T, struct) * diag, slong n,
                               const TEMPLATE(T, ctx_t) ctx)
{
    int parity;
    slong i;
    slong *rows;
    slong *cols;

    rows = _perm_init(mat->r);
    cols = _perm_init(mat->c);

    parity = _perm_randtest(rows, mat->r, state);
    parity ^= _perm_randtest(cols, mat->c, state);

    TEMPLATE(T, mat_zero) (mat, ctx);
    for (i = 0; i < n; i++)
        TEMPLATE(T, mat_entry_set) (mat, rows[i], cols[i], diag + i, ctx);

    _perm_clear(rows);
    _perm_clear(cols);

    return parity;
}


#endif
