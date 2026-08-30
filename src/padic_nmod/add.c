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
padic_nmod_add(padic_nmod_t res, const padic_nmod_t a, const padic_nmod_t b,
               gr_ctx_t ctx)
{
    if (a->u == 0)
    {
        padic_nmod_set(res, b, ctx);
    }

    else if (b->u == 0)
    {
        padic_nmod_set(res, a, ctx);
    }

    else
    {
        if (a->v == b->v)
        {
            res->u = nmod_add(a->u, b->u, PADIC_NMOD_CTX_PN_MOD(ctx));
            res->v = a->v;

            _padic_nmod_canonicalise(res, ctx);
            /* Underflow */
            if (res->v > PADIC_EMAX)
                return GR_UNABLE;
        }

        else if (a->v < b->v)
        {
            ulong f;

            f = (b->v - a->v > PADIC_NMOD_CTX_N(ctx))
                ? PADIC_NMOD_CTX_POW(ctx)[PADIC_NMOD_CTX_N(ctx) - 1]
                : PADIC_NMOD_CTX_POW(ctx)[b->v - a->v - 1];
            res->u = nmod_addmul(a->u, f, b->u, PADIC_NMOD_CTX_PN_MOD(ctx));
            res->v = a->v;
        }

        else  /* a->v > b->v */
        {
            ulong f;

            f = (a->v - b->v > PADIC_NMOD_CTX_N(ctx))
                ? PADIC_NMOD_CTX_POW(ctx)[PADIC_NMOD_CTX_N(ctx) - 1]
                : PADIC_NMOD_CTX_POW(ctx)[a->v - b->v - 1];
            res->u = nmod_addmul(b->u, f, a->u, PADIC_NMOD_CTX_PN_MOD(ctx));
            res->v = b->v;
        }
    }

    return GR_SUCCESS;
}
