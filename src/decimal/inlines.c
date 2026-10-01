/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define DECIMAL_INLINES_C
#include "decimal.h"

void
decball_clear(decball_t res, gr_ctx_t ctx)
{
    decfloat_clear(&res->mid, ctx);
    _decmag_clear(&res->rad, ctx);
}

void
deccfloat_clear(deccfloat_t res, gr_ctx_t ctx)
{
    decfloat_clear(&res->re, ctx);
    decfloat_clear(&res->im, ctx);
}

void
deccball_clear(deccball_t res, gr_ctx_t ctx)
{
    decball_clear(&res->re, ctx);
    decball_clear(&res->im, ctx);
}
