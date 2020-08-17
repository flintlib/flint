/*
    Copyright (C) 2011, 2013 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "padic_mat.h"

void padic_mat_set(padic_mat_t rop, const padic_mat_t op, const padic_ctx_t ctx)
{
    if (op != rop)
    {
        if (padic_mat_val(op) >= padic_mat_prec(rop))
        {
            padic_mat_zero(rop);
        }
        else if (padic_mat_prec(rop) >= padic_mat_prec(op))
        {
            fmpz_mat_set(padic_mat(rop), padic_mat(op));
            padic_mat_val(rop) = padic_mat_val(op);
        }
        else
        {
            fmpz_mat_set(padic_mat(rop), padic_mat(op));
            padic_mat_val(rop) = padic_mat_val(op);

            _padic_mat_reduce(rop, ctx);
        }
    }
}

