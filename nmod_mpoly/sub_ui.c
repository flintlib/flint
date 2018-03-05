/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_sub_ui(nmod_mpoly_t poly1, const nmod_mpoly_t poly2,
                                           ulong c, const nmod_mpoly_ctx_t ctx)
{
    mp_limb_t cr;
    NMOD_RED(cr, c, ctx->ffinfo->mod);
    cr = nmod_neg(cr, ctx->ffinfo->mod);
    nmod_mpoly_add_ui(poly1, poly2, cr, ctx);
}
