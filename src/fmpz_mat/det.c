/*
    Copyright (C) 2010,2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void
fmpz_mat_det(fmpz_t det, const fmpz_mat_t A)
{
    slong dim = A->r;

    if (dim != A->c)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_mat_det). Non-square matrix.\n");
    }

    if (dim < 5)
        fmpz_mat_det_cofactor(det, A);
    else if (dim < 25)
        fmpz_mat_det_bareiss(det, A);
    else if (dim < 60)
        fmpz_mat_det_modular(det, A, 1);
    else
    {
        slong bits = fmpz_mat_max_bits(A);

        if (dim < FLINT_ABS(bits))
            fmpz_mat_det_modular(det, A, 1);
        else
            fmpz_mat_det_modular_accelerated(det, A, 1);
    }
}
