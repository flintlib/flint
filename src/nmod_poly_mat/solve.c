/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "perm.h"

int
nmod_poly_mat_solve(nmod_poly_mat_t X, nmod_poly_t den,
                    const nmod_poly_mat_t A, const nmod_poly_mat_t B)
{
    return nmod_poly_mat_solve_fflu(X, den, A, B);
}
