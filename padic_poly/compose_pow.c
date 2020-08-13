/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "padic_poly.h"

/*
    TODO:  Move this bit of code into "padic".
 */
static void __padic_reduce(fmpz_t u, slong *v, slong N, const padic_ctx_t ctx)
{
    if (!fmpz_is_zero(u))
    {
        if (*v < N)
        {
            int alloc;
            fmpz_t pow;

            alloc = _padic_ctx_pow_ui(pow, N - *v, ctx);
            fmpz_mod(u, u, pow);
            if (alloc)
                fmpz_clear(pow);
        }
        else
        {
            fmpz_zero(u);
            *v = 0;
        }
    }
}

void _padic_poly_compose_pow(fmpz *rop, slong *rval, slong N, 
                             const fmpz *op, slong val, slong len, slong k, 
                             const padic_ctx_t ctx)
{
    if (k == 1)
    {
        if (rop != op)
        {
            _fmpz_vec_set(rop, op, len);
            *rval = val;
        }
    }
    else if (len == 1)
    {
        fmpz_set(rop, op);
        *rval = val;

        __padic_reduce(rop, rval, N, ctx);
    }
    else
    {
        slong i, j, h;

        for (i = len - 1, j = (len - 1) * k ; i >= 0; i--, j -= k)
        {
            fmpz_set(rop + j, op + i);
            if (i)
                for (h = 1; h < k; h++)
                    fmpz_zero(rop + (j - h));
        }
        *rval = val;
    }
}

void padic_poly_compose_pow(padic_poly_t rop, const padic_poly_t op, slong k, 
                            const padic_ctx_t ctx)
{
    const slong len  = op->length;
    const slong lenr = (len - 1) * k + 1;

    if (len == 0)
    {
        padic_poly_zero(rop);
    }
    else
    {
        padic_poly_fit_length(rop, lenr);
        _padic_poly_compose_pow(rop->coeffs, &(rop->val), rop->N, 
                                op->coeffs, op->val, op->length, k, ctx);
        _padic_poly_set_length(rop, lenr);
    }
}

