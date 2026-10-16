/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mpoly.h"

int gr_mpoly_sqr(gr_mpoly_t poly1,
    const gr_mpoly_t poly2,
    gr_mpoly_ctx_t ctx)
{
    /* gr_mpoly_mul detects the repeated operand and dispatches to the
       dedicated squaring routines when the coefficient ring permits it. */
    return gr_mpoly_mul(poly1, poly2, poly2, ctx);
}
