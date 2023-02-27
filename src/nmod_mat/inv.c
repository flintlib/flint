/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_mat.h"


int nmod_mat_inv(nmod_mat_t B, const nmod_mat_t A)
{
    nmod_mat_t I;
    slong i, dim;
    int result;

    dim = A->r;

    switch (dim)
    {
        case 0:
            result = 1;
            break;

        case 1:
            if (nmod_mat_entry(A, 0, 0) == UWORD(0))
            {
                result = 0;
            }
            else
            {
                nmod_mat_entry(B, 0, 0) = 
                    n_invmod(nmod_mat_entry(A, 0, 0), B->mod.n);
                result = 1;
            }
            break;

        default:
            nmod_mat_init(I, dim, dim, B->mod.n);
            for (i = 0; i < dim; i++)
                nmod_mat_entry(I, i, i) = UWORD(1);
            result = nmod_mat_solve(B, A, I);
            nmod_mat_clear(I);
    }

    return result;
}
