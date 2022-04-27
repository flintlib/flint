/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mat.h"

int
fmpq_mat_can_solve(fmpq_mat_t X, const fmpq_mat_t A, const fmpq_mat_t B)
{
    if (fmpq_mat_nrows(A) <= 15)
        return fmpq_mat_can_solve_fraction_free(X, A, B);
    else
        return fmpq_mat_can_solve_multi_mod(X, A, B);
}

