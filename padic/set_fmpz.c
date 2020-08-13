/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic.h"

void padic_set_fmpz(padic_t rop, const fmpz_t op, const padic_ctx_t ctx)
{
    if (!fmpz_is_zero(op))
    {
        padic_val(rop) = fmpz_remove(padic_unit(rop), op, ctx->p);
        _padic_reduce(rop, ctx);
    }
    else
    {
        padic_zero(rop);
    }
}

