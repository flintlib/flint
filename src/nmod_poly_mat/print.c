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
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

void
nmod_poly_mat_print(const nmod_poly_mat_t A, const char * x)
{
    slong i, j;

    flint_printf("<%wd x %wd matrix over Z/nZ[%s]>\n", A->r, A->c, x);

    for (i = 0; i < A->r; i++)
    {
        flint_printf("[");
        for (j = 0; j < A->c; j++)
        {
            /* TODO: pretty */
            nmod_poly_print(nmod_poly_mat_entry(A, i, j));
            if (j + 1 < A->c)
                flint_printf(", ");
        }
        flint_printf("]\n");
    }
    flint_printf("\n");
}
