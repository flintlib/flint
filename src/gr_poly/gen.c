/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

int
gr_poly_gen(gr_poly_t res, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= gr_poly_zero(res, ctx);
    status |= gr_poly_set_coeff_ui(res, 1, 1, ctx);
    return status;
}
