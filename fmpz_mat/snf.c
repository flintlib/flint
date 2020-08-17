/*
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

void
fmpz_mat_snf(fmpz_mat_t S, const fmpz_mat_t A)
{
    fmpz_t det;
    slong m = A->r, n = A->c, b = fmpz_mat_max_bits(A), cutoff = 9;

    if (b <= 2)
        cutoff = 15;
    else if (b <= 4)
        cutoff = 13;
    else if (b <= 8)
        cutoff = 13;
    else if (b <= 16)
        cutoff = 11;
    else if (b <= 32)
        cutoff = 11;
    else if (b <= 64)
        cutoff = 10;

    if (FLINT_MAX(m, n) < cutoff || m != n)
        fmpz_mat_snf_kannan_bachem(S, A);
    else
    {
        fmpz_init(det);
        fmpz_mat_det(det, A);
        if (!fmpz_is_zero(det))
        {
            fmpz_abs(det, det);
            fmpz_mat_snf_iliopoulos(S, A, det);
        }
        else
        {
            fmpz_mat_snf_kannan_bachem(S, A);
        }
        fmpz_clear(det);
    }
}
