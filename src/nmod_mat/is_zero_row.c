/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_mat.h"

int nmod_mat_is_zero_row(const nmod_mat_t mat, slong i)
{
    return _nmod_vec_is_zero(nmod_mat_entry_ptr(mat, i, 0), mat->c);
}
