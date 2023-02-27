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

void
nmod_mat_randfull(nmod_mat_t mat, flint_rand_t state)
{
    slong i;

    for (i = 0; i < mat->r * mat->c; i++)
    {
        mat->entries[i] = FLINT_MAX(1, n_randint(state, mat->mod.n));
    }
}
