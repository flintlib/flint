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
padic_nmod_set_ui(padic_nmod_t res, ulong x, gr_ctx_t ctx)
{
    if (x < PADIC_NMOD_CTX_P(ctx))
    {
        res->man = x;
        res->val = (x == 0) ? PADIC_EMAX : 0;
    }
    else
    {
        ulong p;

        p = PADIC_NMOD_CTX_P(ctx);

        res->val = n_remove2_precomp(&x, p, PADIC_NMOD_CTX_PINV(ctx));
        res->man = nmod_set_ui(x, PADIC_NMOD_CTX_PN_MOD(ctx));
    }

    return GR_SUCCESS;
}
