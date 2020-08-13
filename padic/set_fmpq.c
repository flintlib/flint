/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic.h"

void padic_set_fmpq(padic_t rop, const fmpq_t op, const padic_ctx_t ctx)
{
    if (fmpq_is_zero(op))
    {
        padic_zero(rop);
    }
    else
    {
        fmpq_t t;

        fmpq_init(t);

        padic_val(rop)  = fmpz_remove(fmpq_numref(t), fmpq_numref(op), ctx->p);
        padic_val(rop) -= fmpz_remove(fmpq_denref(t), fmpq_denref(op), ctx->p);

        if (padic_val(rop) >= padic_prec(rop))
        {
            padic_zero(rop);
        }
        else
        {
            _padic_inv(fmpq_denref(t), 
                       fmpq_denref(t), ctx->p, padic_prec(rop) - padic_val(rop));
            fmpz_mul(padic_unit(rop), fmpq_numref(t), fmpq_denref(t));
            _padic_reduce(rop, ctx);
        }
        fmpq_clear(t);
    }
}

