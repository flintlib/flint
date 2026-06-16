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
padic_nmod_set(padic_nmod_t res, const padic_nmod_t x, gr_ctx_t ctx)
{
    res->man = x->man;
    res->val = x->val;

    _padic_nmod_canonicalise(res, ctx);
    /* Underflow */
    if (res->val > PADIC_EMAX)
        return GR_UNABLE;

    return GR_SUCCESS;
}
