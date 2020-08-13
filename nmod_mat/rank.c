/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_mat.h"


slong
nmod_mat_rank(const nmod_mat_t A)
{
    slong m, n, rank;
    slong * perm;
    nmod_mat_t tmp;

    m = A->r;
    n = A->c;

    if (m == 0 || n == 0)
        return 0;

    nmod_mat_init_set(tmp, A);
    perm = flint_malloc(sizeof(slong) * m);

    rank = nmod_mat_lu(perm, tmp, 0);

    flint_free(perm);
    nmod_mat_clear(tmp);
    return rank;
}
