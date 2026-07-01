/*
    Copyright (C) 2026 Rubén Muñoz--Bertrand

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic_nmod.h"

/* Not checked: p is prime */
int
padic_nmod_ctx_init(gr_ctx_t ctx, ulong p, slong n)
{
    slong i;
    ulong hi;

    if (n <= 0 || n >= FLINT_BITS)
        return GR_UNABLE;

    /* Compute p^n and verify that this doesn't overflow a ulong */
    GR_CTX_DATA_AS_PTR(ctx) = flint_malloc(sizeof(_padic_nmod_ctx_struct));
    PADIC_NMOD_CTX_POW(ctx) = (ulong *) flint_malloc(n * sizeof(ulong));

    PADIC_NMOD_CTX_POW(ctx)[0] = p;

    for (i = 1; i < n; i++)
    {
        umul_ppmm(hi, PADIC_NMOD_CTX_POW(ctx)[i],
                  PADIC_NMOD_CTX_POW(ctx)[i - 1], p);

        if (hi != 0)
        {
            padic_nmod_ctx_clear(ctx);
            return GR_UNABLE;
        }
    }

    PADIC_NMOD_CTX_P(ctx) = p;
    PADIC_NMOD_CTX_N(ctx) = n;
    PADIC_NMOD_CTX_PINV1(ctx) = n_binvert(p);
    PADIC_NMOD_CTX_PINV2(ctx) = UWORD_MAX / p;

    nmod_init(&PADIC_NMOD_CTX_P_MOD(ctx), p);
    nmod_init(&PADIC_NMOD_CTX_PN_MOD(ctx), PADIC_NMOD_CTX_POW(ctx)[n - 1]);

    for (i = 0; i < n; i++)
    {
        NMOD_RED(PADIC_NMOD_CTX_POW(ctx)[i], PADIC_NMOD_CTX_POW(ctx)[i],
                 PADIC_NMOD_CTX_PN_MOD(ctx));
    }

    return GR_SUCCESS;
}
