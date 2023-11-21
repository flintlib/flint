/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "gr_mat.h"

/* todo: could reimplement the generic construction here, for
   reduced memory usage in small rings.

   todo: could wrap/reimplement fmpz_mat_is_hadamard */

int
gr_mat_hadamard(gr_mat_t mat, gr_ctx_t ctx)
{
    slong n;
    int status = GR_SUCCESS;
    fmpz_mat_t tmp;

    n = mat->r;

    if (n != mat->c)
        return GR_DOMAIN;

    if (n <= 1)
        return gr_mat_one(mat, ctx);

    if (n > 2 && n % 4 != 0)
        return GR_DOMAIN;

    fmpz_mat_init(tmp, n, n);

    status = fmpz_mat_hadamard(tmp) ? GR_SUCCESS : GR_UNABLE;

    if (status == GR_SUCCESS)
        status = gr_mat_set_fmpz_mat(mat, tmp, ctx);

    fmpz_mat_clear(tmp);

    return status;
}
