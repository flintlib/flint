/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "padic.h"
#include "padic_poly.h"

void _padic_poly_pow(fmpz *rop, slong *rval, slong N,
                     const fmpz *op, slong val, slong len, ulong e,
                     const padic_ctx_t ctx)
{
    fmpz_t pow;
    int alloc;

    *rval = (slong) e * val;

    alloc = _padic_ctx_pow_ui(pow, N - *rval, ctx);

    {
        fmpz * t;
        fmpz_mod_ctx_t mod;

        fmpz_mod_ctx_init(mod, pow);
        t = _fmpz_vec_init(len);

        _fmpz_vec_scalar_mod_fmpz(t, op, len, pow);
        _fmpz_mod_poly_pow(rop, op, len, e, mod);

        fmpz_mod_ctx_clear(mod);
        _fmpz_vec_clear(t, len);
    }

    if (alloc)
        fmpz_clear(pow);
}

void padic_poly_pow(padic_poly_t rop, const padic_poly_t op, ulong e,
                    const padic_ctx_t ctx)
{
    if (e == 0)
    {
        padic_poly_one(rop);
    }
    else if (op->length == 0 || (slong) e * op->val >= rop->N)
    {
        padic_poly_zero(rop);
    }
    else if (e == 1)
    {
        padic_poly_set(rop, op, ctx);
    }
    else
    {
        const slong rlen = (slong) e * (op->length - 1) + 1;
        fmpz *t;

        if (rop == op)
        {
            t = _fmpz_vec_init(rlen);
        }
        else
        {
            padic_poly_fit_length(rop, rlen);
            t = rop->coeffs;
        }

        _padic_poly_pow(t, &(rop->val), rop->N,
                        op->coeffs, op->val, op->length, e, ctx);

        if (rop == op)
        {
            _fmpz_vec_clear(rop->coeffs, rop->alloc);
            rop->coeffs = t;
            rop->alloc  = rlen;
        }
        _padic_poly_set_length(rop, rlen);
        _padic_poly_normalise(rop);
    }
}
