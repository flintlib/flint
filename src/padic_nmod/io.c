/*
    Copyright (C) 2026 Rubén Muñoz--Bertrand

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic_nmod.h"

void
padic_nmod_println(const padic_nmod_t x, gr_ctx_t ctx)
{
    if (x->man == 0)
    {
        flint_printf("0\n");
    }
    else
    {
        flint_printf("%wu*%wu^%wd\n", x->man, PADIC_NMOD_CTX_P(ctx), x->val);
    }
}
