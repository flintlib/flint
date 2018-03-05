/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

void mpoly_ctx_init(mpoly_ctx_t mctx, slong nvars, const ordering_t ord)
{
    mctx->nvars = nvars;
    mctx->ord = ord;
    switch (ord) {
        case ORD_LEX:
            mctx->deg = 0;
            mctx->rev = 0;
            break;
        case ORD_DEGLEX:
            mctx->deg = 1;
            mctx->rev = 0;
            break;
        case ORD_DEGREVLEX:
            mctx->deg = 1;
            mctx->rev = 1;
            break;
        default:
            flint_throw(FLINT_ERROR, "Invalid ordering in mpoly_ctx_init");
    }
    mctx->nfields = mctx->nvars + mctx->deg;
}
