/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mpoly.h"

int gr_mpoly_mul(gr_mpoly_t poly1,
    const gr_mpoly_t poly2,
    const gr_mpoly_t poly3,
    const mpoly_ctx_t mctx, gr_ctx_t cctx)
{
    return gr_mpoly_mul_johnson(poly1, poly2, poly3, mctx, cctx);
}
