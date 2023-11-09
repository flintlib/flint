/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpfr_mat.h"

void
mpfr_mat_swap(mpfr_mat_t mat1, mpfr_mat_t mat2)
{
    FLINT_SWAP(mpfr_mat_struct, *mat1, *mat2);
}
