/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "d_mat.h"

void
d_mat_clear(d_mat_t mat)
{
    if (mat->entries)
    {
        flint_free(mat->entries);   /* Clean up array of entries */
        flint_free(mat->rows);  /* Clean up row array */
    } else if (mat->r != 0)
	flint_free(mat->rows);
}
