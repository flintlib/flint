/*
    Copyright (C) 2011-2012 Fredrik Johansson
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

slong
fmpz_mat_rref(fmpz_mat_t R, fmpz_t den, const fmpz_mat_t A)
{
    if (FLINT_MIN(A->c, A->r) <= 20)
        return fmpz_mat_rref_fflu(R, den, A);
    else if (A->r <= 105 && A->c >= 1.4 * A->r)
        return fmpz_mat_rref_fflu(R, den, A);
    else
        return fmpz_mat_rref_mul(R, den, A);
}

