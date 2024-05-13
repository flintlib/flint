/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2015 Tommy Hofmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define NMOD_MAT_INLINES_C

#include "nmod_mat.h"

void nmod_mat_set_entry(nmod_mat_t mat, slong i, slong j, ulong x)
{
  nmod_mat_entry(mat, i, j) = x;
}
