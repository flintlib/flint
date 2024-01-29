/*
    Copyright (C) 2018 Tommy Hofmann
    Copyright (C) 2020 William Hart
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_mod_mat.h"
#include "gr.h"
#include "gr_mat.h"

int fmpz_mod_mat_can_solve(fmpz_mod_mat_t X, const fmpz_mod_mat_t A,
                           const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    int status;
    _gr_ctx_init_fmpz_mod_from_ref(gr_ctx, ctx);
    status = gr_mat_solve_field((gr_mat_struct *) X, (const gr_mat_struct *) A, (const gr_mat_struct *) B, gr_ctx);
    return (status == GR_SUCCESS) ? 1 : 0;
}
