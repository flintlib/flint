/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz_vec.h"
#include "padic.h"

void padic_ctx_init(padic_ctx_t ctx, const fmpz_t p, slong min, slong max,
                    enum padic_print_mode mode)
{
    if (!(0 <= min && min <= max))
    {
        flint_throw(FLINT_ERROR, "Exception (padic_ctx_init).  Require 0 <= min <= max.");
    }

    fmpz_init_set(ctx->p, p);

    ctx->min = min;
    ctx->max = max;

    ctx->pinv = (!COEFF_IS_MPZ(*p)) ? n_precompute_inverse(fmpz_get_ui(p)) : 0;

    if (max - min > 0)
    {
        slong i, len = max - min;

        ctx->pow = _fmpz_vec_init(len);

        fmpz_pow_ui(ctx->pow, p, ctx->min);
        for (i = 1; i < len; i++)
            fmpz_mul(ctx->pow + i, ctx->pow + (i - 1), p);
    }
    else
    {
        ctx->min = 0;
        ctx->max = 0;
        ctx->pow = NULL;
    }

    ctx->mode = mode;
}

