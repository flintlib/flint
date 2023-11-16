/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_special.h"

/* todo: numerically stable implementation for large re(x) */
int
gr_generic_erfcx(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t;
    int status = GR_SUCCESS;

    GR_TMP_INIT(t, ctx);

    status |= gr_sqr(t, x, ctx);
    status |= gr_exp(t, t, ctx);
    status |= gr_erfc(res, x, ctx);
    status |= gr_mul(res, res, t, ctx);

    GR_TMP_CLEAR(t, ctx);

    if (status != GR_SUCCESS)
        status = GR_UNABLE;

    return status;
}
