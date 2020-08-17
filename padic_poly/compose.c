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

/* Assumes that len1 > 0. */

void _padic_poly_compose(fmpz *rop, slong *rval, slong N, 
                         const fmpz *op1, slong val1, slong len1, 
                         const fmpz *op2, slong val2, slong len2, 
                         const padic_ctx_t ctx)
{
    const slong lenr = (len1 - 1) * (len2 - 1) + 1;

    if (len1 == 1 || len2 == 0)
    {
        fmpz_set(rop, op1);
        *rval = val1;

        __padic_reduce(rop, rval, N, ctx);
    }
    else if (val2 >= 0)
    {
        if (val1 >= N)
        {
            _fmpz_vec_zero(rop, lenr);
            *rval = 0;
        }
        else
        {
            fmpz *vec2;
            fmpz_t f;
            fmpz_t pow;
            int alloc;

            vec2 = _fmpz_vec_init(len2);
            fmpz_init(f);

            fmpz_pow_ui(f, ctx->p, val2);
            _fmpz_vec_scalar_mul_fmpz(vec2, op2, len2, f);

            alloc = _padic_ctx_pow_ui(pow, N - val1, ctx);

            _fmpz_mod_poly_compose(rop, op1, len1, vec2, len2, pow);
            *rval= val1;

            _padic_poly_canonicalise(rop, rval, lenr, ctx->p);

            _fmpz_vec_clear(vec2, len2);
            fmpz_clear(f);
            if (alloc)
                fmpz_clear(pow);
        }
    }
    else  /* val2 < 0 */
    {
        const slong n = len1 - 1;

        if (val1 + n*val2 >= N)
        {
            _fmpz_vec_zero(rop, lenr);
            *rval = 0;
        }
        else
        {
            fmpz_t pow;
            int alloc;
            fmpz *vec1;
            fmpz_t s, t;
            slong i;

            vec1 = _fmpz_vec_init(len1);
            fmpz_init(s);
            fmpz_init(t);

            alloc = _padic_ctx_pow_ui(pow, N - val1 - n*val2, ctx);

            fmpz_pow_ui(s, ctx->p, -val2);
            fmpz_one(t);
            fmpz_set(vec1 + (len1 - 1), op1 + (len1 - 1));
            for (i = len1 - 2; i >= 0; i--)
            {
                fmpz_mul(t, t, s);
                fmpz_mul(vec1 + i, op1 + i, t);
            }

            _fmpz_mod_poly_compose(rop, vec1, len1, op2, len2, pow);
            *rval = val1 + n*val2;

            _padic_poly_canonicalise(rop, rval, lenr, ctx->p);

            _fmpz_vec_clear(vec1, len1);
            fmpz_clear(s);
            fmpz_clear(t);
            if (alloc)
                fmpz_clear(pow);
        }
    }
}

void padic_poly_compose(padic_poly_t rop, 
                        const padic_poly_t op1, const padic_poly_t op2, 
                        const padic_ctx_t ctx)
{
    const slong len1 = op1->length, len2 = op2->length;

    if (len1 == 0)
    {
        padic_poly_zero(rop);
    }
    else if (len1 == 1 || len2 == 0)
    {
        padic_poly_fit_length(rop, 1);
        fmpz_set(rop->coeffs, op1->coeffs);
        rop->val = op1->val;
        _padic_poly_set_length(rop, 1);
        padic_poly_canonicalise(rop, ctx->p);
        padic_poly_reduce(rop, ctx);
    }
    else
    {
        const slong lenr = (len1 - 1) * (len2 - 1) + 1;

        if (rop != op1 && rop != op2)
        {
            padic_poly_fit_length(rop, lenr);
            _padic_poly_compose(rop->coeffs, &(rop->val), rop->N, 
                                op1->coeffs, op1->val, op1->length, 
                                op2->coeffs, op2->val, op2->length, ctx);
            _padic_poly_set_length(rop, lenr);
        }
        else
        {
            fmpz *t = _fmpz_vec_init(lenr);

            _padic_poly_compose(t, &(rop->val), rop->N, 
                                op1->coeffs, op1->val, op1->length, 
                                op2->coeffs, op2->val, op2->length, ctx);
            _fmpz_vec_clear(rop->coeffs, rop->alloc);
            rop->coeffs = t;
            rop->alloc  = lenr;
            rop->length = lenr;
        }
        _padic_poly_normalise(rop);
    }
}

