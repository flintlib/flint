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
ca_mat_lu(slong * rank, slong * P, ca_mat_t LU, const ca_mat_t A, int rank_check, ca_ctx_t ctx)
{
    return ca_mat_lu_recursive(rank, P, LU, A, rank_check, ctx);
}
