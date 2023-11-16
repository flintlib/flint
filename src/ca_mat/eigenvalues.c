/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

int
ca_mat_eigenvalues(ca_vec_t lambda, ulong * exp, const ca_mat_t mat, ca_ctx_t ctx)
{
    int success;
    ca_poly_t cp;
    ca_poly_init(cp, ctx);
    ca_mat_charpoly(cp, mat, ctx);
    success = ca_poly_roots(lambda, exp, cp, ctx);
    ca_poly_clear(cp, ctx);
    return success;
}
