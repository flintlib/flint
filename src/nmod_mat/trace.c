/*
    Copyright (C) 2011 Fredrik Johansson

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

mp_limb_t
nmod_mat_trace(const nmod_mat_t mat)
{
    mp_limb_t t;
    slong i, n = nmod_mat_nrows(mat);

    if (n == 0)
        return 0;

    t = nmod_mat_entry(mat, 0, 0);

    for (i = 1; i < n; i++)
        t = nmod_add(t, nmod_mat_entry(mat, i, i), mat->mod);

    return t;
}
