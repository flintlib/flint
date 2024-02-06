/*
    Copyright (C) 2018 Daniel Schultz
    Copyright (C) 2020, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "fmpz_mpoly.h"
#include "qqbar.h"
#include "gr.h"
#include "gr_generic.h"

int
qqbar_evaluate_fmpz_mpoly(qqbar_t res, const fmpz_mpoly_t f, qqbar_srcptr x, slong deg_limit, slong bits_limit, const fmpz_mpoly_ctx_t ctx)
{
    gr_ctx_t QQbar;
    gr_ctx_init_complex_qqbar(QQbar);
    _gr_ctx_qqbar_set_limits(QQbar, deg_limit, bits_limit);

    if (f->length <= 1)
        return (gr_fmpz_mpoly_evaluate_iter(res, f, x, ctx, QQbar) == GR_SUCCESS);
    else
        return (gr_fmpz_mpoly_evaluate_horner(res, f, x, ctx, QQbar) == GR_SUCCESS);
}
