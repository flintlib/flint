/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_factor_insert(ca_factor_t fac, const ca_t base, const ca_t exp, ca_ctx_t ctx)
{
    slong i;

    for (i = 0; i < fac->length; i++)
    {
        if (ca_equal_repr(fac->base + i, base, ctx))
        {
            ca_add(fac->exp, fac->exp, exp, ctx);
            return;
        }
    }

    if (fac->length == fac->alloc)
    {
        slong new_alloc;

        new_alloc = FLINT_MAX(1, 2 * fac->alloc);

        fac->base = (ca_ptr) flint_realloc(fac->base, sizeof(ca_struct) * new_alloc);
        fac->exp = (ca_ptr) flint_realloc(fac->exp, sizeof(ca_struct) * new_alloc);

        for (i = fac->alloc; i < new_alloc; i++)
        {
            ca_init(fac->base + i, ctx);
            ca_init(fac->exp + i, ctx);
        }

        fac->alloc = new_alloc;
    }

    ca_set(fac->base + fac->length, base, ctx);
    ca_set(fac->exp + fac->length, exp, ctx);
    fac->length++;
}

