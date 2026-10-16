/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "padic.h"
#include "padic_poly.h"

void _padic_poly_mul(fmpz *rop, slong *rval, slong N,
                     const fmpz *op1, slong val1, slong len1,
                     const fmpz *op2, slong val2, slong len2,
                     const padic_ctx_t ctx)
{
    fmpz_t pow;
    int alloc;

    *rval = val1 + val2;

    alloc = _padic_ctx_pow_ui(pow, N - *rval, ctx);

    _fmpz_poly_mul(rop, op1, len1, op2, len2);
    _fmpz_vec_scalar_mod_fmpz(rop, rop, len1 + len2 - 1, pow);

    if (alloc)
        fmpz_clear(pow);
}

void padic_poly_mul(padic_poly_t res,
                    const padic_poly_t poly1, const padic_poly_t poly2,
                    const padic_ctx_t ctx)
{
    const slong lenG = poly1->length;
    const slong lenH = poly2->length;
    const slong lenF = lenG + lenH - 1;

    if (lenG == 0 || lenH == 0 || poly1->val + poly2->val >= res->N)
    {
        padic_poly_zero(res);
    }
    else
    {
        fmpz *t;

        if (res == poly1 || res == poly2)
        {
            t = _fmpz_vec_init(lenF);
        }
        else
        {
            padic_poly_fit_length(res, lenF);
            t = res->coeffs;
        }

        if (lenG >= lenH)
            _padic_poly_mul(t, &(res->val), res->N, poly1->coeffs, poly1->val, lenG,
                                          poly2->coeffs, poly2->val, lenH, ctx);
        else
            _padic_poly_mul(t, &(res->val), res->N, poly2->coeffs, poly2->val, lenH,
                                          poly1->coeffs, poly1->val, lenG, ctx);

        if (res == poly1 || res == poly2)
        {
            _fmpz_vec_clear(res->coeffs, res->alloc);
            res->coeffs = t;
            res->alloc  = lenF;
        }

        _padic_poly_set_length(res, lenF);
        _padic_poly_normalise(res);
    }
}
