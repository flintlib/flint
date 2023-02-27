/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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

void padic_poly_mul(padic_poly_t f, 
                    const padic_poly_t g, const padic_poly_t h, 
                    const padic_ctx_t ctx)
{
    const slong lenG = g->length;
    const slong lenH = h->length;
    const slong lenF = lenG + lenH - 1;

    if (lenG == 0 || lenH == 0 || g->val + h->val >= f->N)
    {
        padic_poly_zero(f);
    }
    else
    {
        fmpz *t;

        if (f == g || f == h)
        {
            t = _fmpz_vec_init(lenF);
        }
        else
        {
            padic_poly_fit_length(f, lenF);
            t = f->coeffs;
        }

        if (lenG >= lenH)
            _padic_poly_mul(t, &(f->val), f->N, g->coeffs, g->val, lenG, 
                                          h->coeffs, h->val, lenH, ctx);
        else
            _padic_poly_mul(t, &(f->val), f->N, h->coeffs, h->val, lenH, 
                                          g->coeffs, g->val, lenG, ctx);

        if (f == g || f == h)
        {
            _fmpz_vec_clear(f->coeffs, f->alloc);
            f->coeffs = t;
            f->alloc  = lenF;
        }

        _padic_poly_set_length(f, lenF);
        _padic_poly_normalise(f);
    }
}

