/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

void
fmpz_poly_mat_clear(fmpz_poly_mat_t A)
{
    if (A->entries)
    {
        slong i;

        for (i = 0; i < A->r * A->c; i++)
            fmpz_poly_clear(A->entries + i);

        flint_free(A->entries);
        flint_free(A->rows);
    } else if (A->r != 0)
        flint_free(A->rows);
}
