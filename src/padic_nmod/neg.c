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
padic_nmod_neg(padic_nmod_t res, const padic_nmod_t a, gr_ctx_t ctx)
{
    if (a->man == 0)
    {
        padic_nmod_zero(res, ctx);
    }
    else
    {
        res->val = a->val;
        res->man = nmod_neg(a->man, PADIC_NMOD_CTX_PN_MOD(ctx));
    }

    return GR_SUCCESS;
}
