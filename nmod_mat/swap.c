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

void
nmod_mat_swap(nmod_mat_t mat1, nmod_mat_t mat2)
{
    nmod_mat_t temp;
    *temp = *mat1;
    *mat1 = *mat2;
    *mat2 = *temp;
}
