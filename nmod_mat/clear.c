/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"

void
nmod_mat_clear(nmod_mat_t mat)
{
    if (mat->entries)
        flint_free(mat->entries);
    if (mat->rows)
        flint_free(mat->rows);
    memset(mat, 0, sizeof(*mat));
}
