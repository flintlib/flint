/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mat.h"

slong fmpz_mod_mat_rank(const fmpz_mod_mat_t A)
{
    slong m, n, rank;
    slong *perm;
    fmpz_mod_mat_t tmp;

    m = A->mat->r;
    n = A->mat->c;

    if (m == 0 || n == 0)
        return 0;

    fmpz_mod_mat_init_set(tmp, A);
    perm = flint_malloc(sizeof(slong) * m);

    rank = fmpz_mod_mat_lu(perm, tmp, 0);

    flint_free(perm);
    fmpz_mod_mat_clear(tmp);
    return rank;
}

