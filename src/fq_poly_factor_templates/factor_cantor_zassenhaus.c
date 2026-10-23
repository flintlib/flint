/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "gr.h"
#include "gr_poly.h"
#include "templates.h"

void
TEMPLATE(T, poly_factor_cantor_zassenhaus) (TEMPLATE(T, poly_factor_t) res,
    const TEMPLATE(T, poly_t) f, const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_factor_gr_algorithm) (res, NULL, f,
        GR_POLY_FACTOR_ALGORITHM_CANTOR_ZASSENHAUS, ctx);
}

#endif
