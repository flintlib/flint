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
nmod_mat_neg(nmod_mat_t B, const nmod_mat_t A)
{
    slong i;

    if (A->c == 0)
        return;

    for (i = 0; i < A->r; i++)
        _nmod_vec_neg(B->rows[i], A->rows[i], A->c, A->mod);
}
