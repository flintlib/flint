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
ca_mat_conj_transpose(ca_mat_t mat1, const ca_mat_t mat2, ca_ctx_t ctx)
{
    ca_mat_transpose(mat1, mat2, ctx);
    ca_mat_conj(mat1, mat1, ctx);
}
