/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_poly.h"
#include "fmpq_mat.h"
#include "qqbar.h"

void
qqbar_eigenvalues_fmpq_mat(qqbar_ptr res, const fmpq_mat_t mat, int flags)
{
    fmpq_poly_t t;
    fmpq_poly_init(t);
    fmpq_mat_charpoly(t, mat);
    qqbar_roots_fmpq_poly(res, t, flags);
    fmpq_poly_clear(t);
}
