/*
    Copyright (C) 2026 Rubén Muñoz--Bertrand

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "long_extras.h"
#include "padic_nmod.h"

#define PADIC_NMOD_RANDTEST_TRIES 10

int
padic_nmod_randtest(padic_nmod_t rop, flint_rand_t state, gr_ctx_t ctx)
{
    rop->man = _n_randint(state, PADIC_NMOD_CTX_PN_MOD(ctx).n);

    if (!rop->man)
    {
        rop->val = 0;
    }
    else
    {
        rop->val = z_randint(state, PADIC_EMAX + 1);
        _padic_nmod_canonicalise(rop, ctx);
    }

    return GR_SUCCESS;
}

int
padic_nmod_randtest_not_zero(padic_nmod_t rop, flint_rand_t state,
                             gr_ctx_t ctx)
{
    slong i = 1;

    do
    {
        padic_nmod_randtest(rop, state, ctx);
        i++;
    }
    while (padic_nmod_is_zero(rop, ctx) == T_TRUE
           && i < PADIC_NMOD_RANDTEST_TRIES);

    if (padic_nmod_is_zero(rop, ctx) == T_TRUE)
    {
        rop->man = 1;
        rop->val = PADIC_EMAX;
    }

    return GR_SUCCESS;
}
