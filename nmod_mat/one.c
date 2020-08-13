/*
    Copyright (C) 2014 Ashish Kedia

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

/*Function to convert a square matrix to an identity matrix
    The matrix is assumed to be a square one*/
void
nmod_mat_one(nmod_mat_t mat)
{
    slong i,j;
    for(i = 0; i < mat->r; i++)
        for(j = 0; j < mat->c; j++)
            if(i==j) nmod_mat_entry(mat, i, j) = 1;
            else nmod_mat_entry(mat, i, j) = 0;
}
