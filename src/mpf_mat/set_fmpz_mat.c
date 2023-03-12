/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "mpf_mat.h"
#include "fmpz_vec.h"

void
mpf_mat_set_fmpz_mat(mpf_mat_t B, const fmpz_mat_t A)
{
    slong i;

    if (A->c != 0)
        for (i = 0; i < A->r; i++)
            _mpf_vec_set_fmpz_vec(B->rows[i], A->rows[i], A->c);
}
