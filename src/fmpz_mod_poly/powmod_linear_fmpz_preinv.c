/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"


static void _fmpz_mod_poly_powmod_linear_fmpz_preinv(
    fmpz * res,
    const fmpz_t a,
    const fmpz_t e,
    const fmpz * f, slong lenf,
    const fmpz* finv, slong lenfinv,
    const fmpz_mod_ctx_t ctx)
{
    fmpz * T, * Q;
    slong lenT = 2 * lenf - 3;
    slong lenQ = lenT - lenf + 1;
    slong i, j;
    fmpz_t t, lcinv;

    FLINT_ASSERT(lenf > 1);
    FLINT_ASSERT(!fmpz_is_zero(e));

    fmpz_init(t);

    if (lenf <= 2)
    {
        fmpz_mod_mul(t, f + 0, finv + 0, ctx);
        fmpz_mod_sub(t, a, t, ctx);
        fmpz_mod_pow_fmpz(res + 0, t, e, ctx);
        fmpz_clear(t);
        return;
    }

    fmpz_init(lcinv);
    T = _fmpz_vec_init(lenT + lenQ);
    Q = T + lenT;

    _fmpz_vec_zero(res, lenf - 1);
    fmpz_set(res + 0, a);
    fmpz_one(res + 1);

    for (i = fmpz_sizeinbase(e, 2) - 2; i >= 0; i--)
    {
        _fmpz_mod_poly_sqr(T, res, lenf - 1, fmpz_mod_ctx_modulus(ctx));
        _fmpz_mod_poly_divrem_newton_n_preinv(Q, res, T, 2 * lenf - 3, f, lenf,
                                     finv, lenfinv, fmpz_mod_ctx_modulus(ctx));
        if (fmpz_tstbit(e, i))
        {
            j = lenf - 1;
            fmpz_mod_mul(lcinv, finv + 0, res + j - 1, ctx);
            fmpz_mod_neg(lcinv, lcinv, ctx);
            for (j--; j > 0; j--)
            {
                fmpz_mul(t, a, res + j);
                fmpz_addmul(t, lcinv, f + j);
                fmpz_add(t, t, res + j - 1);
                fmpz_mod(res + j, t, fmpz_mod_ctx_modulus(ctx));
            }
            fmpz_mul(t, a, res + j);
            fmpz_addmul(t, lcinv, f + j);
            fmpz_mod(res + j, t, fmpz_mod_ctx_modulus(ctx));
        }
    }

    fmpz_clear(lcinv);
    fmpz_clear(t);
    _fmpz_vec_clear(T, lenT + lenQ);
}


/* res = (x+a)^e mod f */
void fmpz_mod_poly_powmod_linear_fmpz_preinv(
    fmpz_mod_poly_t res,
    const fmpz_t a,
    const fmpz_t e,
    const fmpz_mod_poly_t f,
    const fmpz_mod_poly_t finv,
    const fmpz_mod_ctx_t ctx)
{
    slong lenf = f->length;
    slong trunc = lenf - 1;
    int sgn = fmpz_sgn(e);
    fmpz_mod_poly_t tmp;

    if (lenf < 2)
    {
        fmpz_mod_poly_zero(res, ctx);
        return;
    }

    if (sgn < 0)
    {
        flint_throw(FLINT_ERROR, "fmpz_mod_poly_powmod_linear_fmpz_preinv: "
                                               "negative exp not implemented");
    }

    if (sgn == 0)
    {
        fmpz_mod_poly_one(res, ctx);
        return;
    }

    if (res == f || res == finv)
    {
        fmpz_mod_poly_init2(tmp, trunc, ctx);
        _fmpz_mod_poly_powmod_linear_fmpz_preinv(tmp->coeffs, a, e,
                             f->coeffs, lenf, finv->coeffs, finv->length, ctx);
        fmpz_mod_poly_swap(res, tmp, ctx);
        fmpz_mod_poly_clear(tmp, ctx);
    }
    else
    {
        fmpz_mod_poly_fit_length(res, trunc, ctx);
        _fmpz_mod_poly_powmod_linear_fmpz_preinv(res->coeffs, a, e,
                             f->coeffs, lenf, finv->coeffs, finv->length, ctx);
    }

    _fmpz_mod_poly_set_length(res, trunc);
    _fmpz_mod_poly_normalise(res);
}

