/*
    Copyright (C) 2011-2012, 2025 Fredrik Johansson
    Copyright (C) 2014 Alex J. Best

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

slong
fmpz_mat_rref(fmpz_mat_t R, fmpz_t den, const fmpz_mat_t A)
{
    slong r = A->r;
    slong c = A->c;

    if (r <= 3 || c <= 2 || (r <= 20 && c > r) || (r > 20 && r <= 100 && c > r + (r - 20) / 80.0 * r))
        return fmpz_mat_rref_fflu(R, den, A);
    else
        return fmpz_mat_rref_mul(R, den, A);
}
