/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"
#include "gr.h"
#include "gr_mat.h"

int
fmpz_mat_is_one(const fmpz_mat_t mat)
{
    gr_ctx_t ctx;
    gr_ctx_init_fmpz(ctx);
    return gr_mat_is_one((const gr_mat_struct *) mat, ctx) == T_TRUE;
}
