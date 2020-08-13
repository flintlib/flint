/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "padic_mat.h"

void padic_mat_set_entry_padic(padic_mat_t rop, slong i, slong j, 
                               const padic_t op, const padic_ctx_t ctx)
{
    if (padic_is_zero(op))
    {
        fmpz_zero(padic_mat_entry(rop, i, j));
        _padic_mat_canonicalise(rop, ctx);
        return;
    }

    if (padic_mat_val(rop) == padic_val(op))
    {
        fmpz_set(padic_mat_entry(rop, i, j), padic_unit(op));
    }
    else if (padic_mat_val(rop) < padic_val(op))
    {
        fmpz_t t;

        fmpz_init(t);
        fmpz_pow_ui(t, ctx->p, padic_val(op) - padic_mat_val(rop));
        fmpz_mul(padic_mat_entry(rop, i, j), padic_unit(op), t);
        fmpz_clear(t);

        _padic_mat_canonicalise(rop, ctx);
    }
    else  /* rop->val > op->val */
    {
        fmpz_t pow;

        fmpz_init(pow);
        fmpz_pow_ui(pow, ctx->p, padic_mat_val(rop) - padic_val(op));
        _fmpz_vec_scalar_mul_fmpz(padic_mat(rop)->entries, 
                                  padic_mat(rop)->entries, 
                                  padic_mat(rop)->r * padic_mat(rop)->c, pow);
        fmpz_clear(pow);

        fmpz_set(padic_mat_entry(rop, i, j), padic_unit(op));
        padic_mat_val(rop) = padic_val(op);
    }
}

