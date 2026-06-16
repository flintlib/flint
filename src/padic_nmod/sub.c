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
padic_nmod_sub(padic_nmod_t res, const padic_nmod_t a, const padic_nmod_t b,
               gr_ctx_t ctx)
{
    if (a->man == 0)
    {
        padic_nmod_neg(res, b, ctx);
    }
    else if (b->man == 0)
    {
        padic_nmod_set(res, a, ctx);
    }
    else
    {
        if (a->val == b->val)
        {
            res->man = nmod_sub(a->man, b->man, PADIC_NMOD_CTX_PN_MOD(ctx));
            res->val = a->val;

            _padic_nmod_canonicalise(res, ctx);
            /* Underflow */
            if (res->val > PADIC_EMAX)
                return GR_UNABLE;
        }
        else if (a->val < b->val)
        {
            nmod_t pn = PADIC_NMOD_CTX_PN_MOD(ctx);
            ulong f, n;

            f = (b->val - a->val > PADIC_NMOD_CTX_N(ctx))
                ? PADIC_NMOD_CTX_POW(ctx)[PADIC_NMOD_CTX_N(ctx) - 1]
                : PADIC_NMOD_CTX_POW(ctx)[b->val - a->val - 1];
            n = nmod_neg(b->man, pn);
            res->man = nmod_addmul(a->man, f, n, pn);
            res->val = a->val;
        }
        else  /* a->val > b->val */
        {
            nmod_t pn = PADIC_NMOD_CTX_PN_MOD(ctx);
            ulong f;

            f = (a->val - b->val > PADIC_NMOD_CTX_N(ctx))
                ? PADIC_NMOD_CTX_POW(ctx)[PADIC_NMOD_CTX_N(ctx) - 1]
                : PADIC_NMOD_CTX_POW(ctx)[a->val - b->val - 1];
            res->man = nmod_neg(b->man, pn);
            res->man = nmod_addmul(res->man, f, a->man, pn);
            res->val = b->val;
        }
    }

    return GR_SUCCESS;
}
