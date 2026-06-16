/*
    Copyright (C) 2026 Rubén Muñoz--Bertrand

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic_nmod.h"

int
padic_nmod_div(padic_nmod_t res, const padic_nmod_t a, const padic_nmod_t b,
               gr_ctx_t ctx)
{
    if (b->man == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (padic_nmod_div).  b is zero.\n");
        return GR_UNABLE;
    }

    if (a->man == 0)
        return padic_nmod_zero(res, ctx);

    res->man = nmod_div(a->man, b->man, PADIC_NMOD_CTX_PN_MOD(ctx));
    res->val = a->val - b->val;

    /* Overflow or underflow */
    if (res->val < PADIC_EMIN || res->val > PADIC_EMAX)
        return GR_UNABLE;

    return GR_SUCCESS;
}
