/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "mpf_mat.h"

void
mpf_mat_clear(mpf_mat_t mat)
{
    if (mat->entries)
    {
        slong i;
        for (i = 0; i < mat->r * mat->c; i++)
            mpf_clear(mat->entries + i);    /* Clear all coefficients */
        flint_free(mat->entries);   /* Clean up array of entries */
        flint_free(mat->rows);  /* Clean up row array */
    }
}
