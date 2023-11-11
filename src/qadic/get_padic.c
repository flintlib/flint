/*
    Copyright (C) 2011, 2012, 2013 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "padic.h"
#include "qadic.h"

int qadic_get_padic(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
{
    if (op->length > 0)
    {
        if (_fmpz_vec_is_zero(op->coeffs + 1, op->length - 1))
        {
            fmpz_set(padic_unit(rop), op->coeffs + 0);
            padic_val(rop) = op->val;
            _padic_canonicalise(rop, &ctx->pctx);
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        padic_zero(rop);
        return 1;
    }
}
