/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void
fmpz_mat_get_mpf_mat(mpf_mat_t B, const fmpz_mat_t A)
{
    slong i;

    if (A->c != 0)
        for (i = 0; i < A->r; i++)
            _fmpz_vec_get_mpf_vec(B->rows[i], A->rows[i], A->c);
}
