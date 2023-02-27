/*
    Copyright (C) 2015 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "ulong_extras.h"

slong
nmod_mat_howell_form(nmod_mat_t A)
{
    slong i, j, n;
    slong k;

    n = A->r;
    k = n;

    if (nmod_mat_is_empty(A))
        return 0;

    nmod_mat_strong_echelon_form(A);

    for (i = 0; i < n; i++)
    {
        if (nmod_mat_is_zero_row(A, i))
        {
            k--;
            for (j = i + 1; j < n; j++)
            {
                if (!nmod_mat_is_zero_row(A, j))
                {
                    nmod_mat_swap_rows(A, NULL, i, j);
                    j = n;
                    k++;
                }
            }
        }
    }
    return k;
}

