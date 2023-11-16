/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "fmpz_poly.h"
#include "qqbar.h"

void
qqbar_eigenvalues_fmpz_mat(qqbar_ptr res, const fmpz_mat_t mat, int flags)
{
    fmpz_poly_t t;
    fmpz_poly_init(t);
    fmpz_mat_charpoly(t, mat);
    qqbar_roots_fmpz_poly(res, t, flags);
    fmpz_poly_clear(t);
}
