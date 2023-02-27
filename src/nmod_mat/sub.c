/*
    Copyright (C) 2010 Fredrik Johansson

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
nmod_mat_sub(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    slong i;

    if (C->c == 0)
        return;

    for (i = 0; i < C->r; i++)
    {
        _nmod_vec_sub(C->rows[i], A->rows[i], B->rows[i], C->c, C->mod);
    }
}
