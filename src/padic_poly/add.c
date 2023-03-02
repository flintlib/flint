/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "padic_poly.h"

void _padic_poly_add(fmpz *rop, slong *val, slong N, 
                     const fmpz *op1, slong val1, slong len1, slong N1, 
                     const fmpz *op2, slong val2, slong len2, slong N2, 
                     const padic_ctx_t ctx)
{
    const slong len = FLINT_MAX(len1, len2);

    *val = FLINT_MIN(val1, val2);

    if (val1 == val2)
    {
        _fmpz_poly_add(rop, op1, len1, op2, len2);
        _padic_poly_canonicalise(rop, val, len, ctx->p);
    }
    else  /* => (op1 != op2) */
    {
        fmpz_t x;

        fmpz_init(x);
        if (val1 < val2)  /* F := p^g (G + p^{h-g} H) */
        {
            fmpz_pow_ui(x, ctx->p, val2 - val1);

            if (rop == op1)
            {
                _fmpz_vec_zero(rop + len1, len2 - len1);
                _fmpz_vec_scalar_addmul_fmpz(rop, op2, len2, x);
            }
            else
            {
                _fmpz_vec_scalar_mul_fmpz(rop, op2, len2, x);
                _fmpz_poly_add(rop, op1, len1, rop, len2);
            }
        }
        else  /* F := p^h (p^{g-h} G + H) */
        {
            fmpz_pow_ui(x, ctx->p, val1 - val2);

            if (rop == op2)
            {
                _fmpz_vec_zero(rop + len2, len1 - len2);
                _fmpz_vec_scalar_addmul_fmpz(rop, op1, len1, x);
            }
            else
            {
                _fmpz_vec_scalar_mul_fmpz(rop, op1, len1, x);
                _fmpz_poly_add(rop, rop, len1, op2, len2);
            }
        }
        fmpz_clear(x);
    }

    /* Reduce */
    if (N - *val > 0)
    {
        fmpz_t pow;
        int alloc;

        alloc = _padic_ctx_pow_ui(pow, N - *val, ctx);

        _fmpz_vec_scalar_mod_fmpz(rop, rop, len, pow);

        if (_fmpz_vec_is_zero(rop, len))
            *val = 0;

        if (alloc)
            fmpz_clear(pow);
    }
    else
    {
        _fmpz_vec_zero(rop, len);
        *val = 0;
    }
}

void padic_poly_add(padic_poly_t f, 
                    const padic_poly_t g, const padic_poly_t h, 
                    const padic_ctx_t ctx)
{
    const slong lenG = g->length;
    const slong lenH = h->length;
    const slong lenF = FLINT_MAX(lenG, lenH);

    if (lenG == 0)
    {
        padic_poly_set(f, h, ctx);
        return;
    }
    if (lenH == 0)
    {
        padic_poly_set(f, g, ctx);
        return;
    }
    if ((lenG == 0 && lenH == 0) || (FLINT_MIN(g->val, h->val) >= f->N))
    {
        padic_poly_zero(f);
        return;
    }

    padic_poly_fit_length(f, lenF);

    _padic_poly_add(f->coeffs, &(f->val), f->N, 
                    g->coeffs, g->val, lenG, g->N, 
                    h->coeffs, h->val, lenH, h->N, ctx);

    _padic_poly_set_length(f, lenF);
    _padic_poly_normalise(f);
}

