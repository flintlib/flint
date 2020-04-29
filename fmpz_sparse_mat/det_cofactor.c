/*
    Copyright (C) 2010,2011,2018 Fredrik Johansson

    This file is part of FLINT.
    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_mat.h"

void
fmpz_sparse_mat_det_cofactor(fmpz_t det, const fmpz_sparse_mat_t M)
{
    fmpz_mat_t dM;
    if (M->r != M->c) {fmpz_zero(det); return;}
    if (M->r > 4)
    {
        flint_printf("Exception (fmpz_sparse_mat_det_cofactor). dim > 4 not implemented.");
        flint_abort();
    }
    fmpz_mat_init(dM, M->r, M->c);
    fmpz_sparse_mat_to_dense(dM, M);
    fmpz_mat_det_cofactor(det, dM);
    fmpz_mat_clear(dM);
}
