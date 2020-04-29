/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_mat.h"
#include "ulong_extras.h"

void
nmod_sparse_mat_print_pretty(const nmod_sparse_mat_t M)
{
    slong i;
    char row_fmt[FLINT_BITS + 5];
    flint_sprintf(row_fmt, "%%%dwd: ", n_sizeinbase(M->r, 10));

    flint_printf("<%wd x %wd sparse integer matrix mod %w>\n", 
                M->r, M->c, M->mod.n);

    for (i = 0; i < M->r; i++)
    {
        flint_printf(row_fmt, i);
        nmod_sparse_vec_print_pretty(&M->rows[i], M->c_off, M->c, M->mod);
    }
}

