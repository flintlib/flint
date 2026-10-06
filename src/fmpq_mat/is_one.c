/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpq_mat.h"
#include "gr.h"
#include "gr_mat.h"

int
fmpq_mat_is_one(const fmpq_mat_t mat)
{
    gr_ctx_t ctx;
    gr_ctx_init_fmpq(ctx);
    return gr_mat_is_one((const gr_mat_struct *) mat, ctx) == T_TRUE;
}
