/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

void
ca_mat_det_berkowitz(ca_t res, const ca_mat_t A, ca_ctx_t ctx)
{
    ca_ptr t;
    t = _ca_vec_init(ca_mat_nrows(A) + 1, ctx);

    _ca_mat_charpoly(t, A, ctx);
    ca_swap(res, t, ctx);
    if (ca_mat_nrows(A) % 2)
        ca_neg(res, res, ctx);

    _ca_vec_clear(t, ca_mat_nrows(A) + 1, ctx);
}
