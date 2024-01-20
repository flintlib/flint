/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2015 Vladimir Glazachev
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_poly.h"

void fmpz_mod_poly_scalar_addmul_fmpz(fmpz_mod_poly_t A,
             const fmpz_mod_poly_t B, const fmpz_t x, const fmpz_mod_ctx_t ctx)
{
    slong len = FLINT_MAX(A->length, B->length);

    if (fmpz_is_zero(x) || B->length < 1)
        return;

    fmpz_mod_poly_fit_length(A, B->length, ctx);

    if (B->length > A->length)
        _fmpz_vec_zero(A->coeffs + A->length, B->length - A->length);

    _fmpz_vec_scalar_addmul_fmpz(A->coeffs, B->coeffs, B->length, x);
    _fmpz_mod_vec_set_fmpz_vec(A->coeffs, A->coeffs, len, ctx);

    _fmpz_mod_poly_set_length(A, len);
    _fmpz_mod_poly_normalise(A);
}

void _fmpz_mod_poly_scalar_div_fmpz(fmpz *res, const fmpz *poly, slong len,
                                    const fmpz_t x, const fmpz_mod_ctx_t ctx)
{
    fmpz_t g, xinv;

    fmpz_init(g);
    fmpz_init(xinv);

    if (fmpz_sgn(x) < 0 || fmpz_cmp(x, fmpz_mod_ctx_modulus(ctx)) >= 0)
    {
       fmpz_mod(xinv, x, fmpz_mod_ctx_modulus(ctx));
       fmpz_gcdinv(g, xinv, xinv, fmpz_mod_ctx_modulus(ctx));
    } else
       fmpz_gcdinv(g, xinv, x, fmpz_mod_ctx_modulus(ctx));

    if (!fmpz_is_one(g))
    {
        flint_throw(FLINT_ERROR, "Exception (_fmpz_mod_poly_scalar_div_fmpz). Impossible inverse.\n");
    }

    _fmpz_mod_vec_scalar_mul_fmpz_mod(res, poly, len, xinv, ctx);

    fmpz_clear(xinv);
    fmpz_clear(g);
}

void fmpz_mod_poly_scalar_div_fmpz(fmpz_mod_poly_t res,
          const fmpz_mod_poly_t poly, const fmpz_t x, const fmpz_mod_ctx_t ctx)
{

    if (fmpz_is_zero(x))
    {
        if (fmpz_is_one(fmpz_mod_ctx_modulus(ctx)))
        {
            fmpz_mod_poly_set(res, poly, ctx);
            return;
        }
        else
        {
            flint_throw(FLINT_ERROR, "Exception (fmpz_mod_poly_scalar_div_fmpz). Division by zero.\n");
        }
    }

    fmpz_mod_poly_fit_length(res, poly->length, ctx);
    _fmpz_mod_poly_scalar_div_fmpz(res->coeffs, poly->coeffs, poly->length, x, ctx);

    _fmpz_mod_poly_set_length(res, poly->length);
    _fmpz_mod_poly_normalise(res);
}

void _fmpz_mod_poly_scalar_mul_fmpz(fmpz *res, const fmpz *poly, slong len,
                                    const fmpz_t x, const fmpz_mod_ctx_t ctx)
{
    if (fmpz_sgn(x) >= 0 && fmpz_cmp(x, fmpz_mod_ctx_modulus(ctx)) < 0)
    {
        _fmpz_mod_vec_scalar_mul_fmpz_mod(res, poly, len, x, ctx);
    }
    else
    {
        fmpz_t t;
        fmpz_init(t);

        /* todo: handle small negative coefficients specially */
        fmpz_mod_set_fmpz(t, x, ctx);

        _fmpz_mod_vec_scalar_mul_fmpz_mod(res, poly, len, t, ctx);
        fmpz_clear(t);
    }
}

void fmpz_mod_poly_scalar_mul_fmpz(fmpz_mod_poly_t res,
          const fmpz_mod_poly_t poly, const fmpz_t x, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_fit_length(res, poly->length, ctx);
    _fmpz_mod_poly_scalar_mul_fmpz(res->coeffs, poly->coeffs, poly->length, x, ctx);
    _fmpz_mod_poly_set_length(res, poly->length);
    _fmpz_mod_poly_normalise(res);
}

void _fmpz_mod_poly_scalar_mul_ui(fmpz *res, const fmpz *poly, slong len,
                                    ulong x, const fmpz_mod_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_mod_set_ui(t, x, ctx);
    _fmpz_mod_vec_scalar_mul_fmpz_mod(res, poly, len, t, ctx);
    fmpz_clear(t);
}

void fmpz_mod_poly_scalar_mul_ui(fmpz_mod_poly_t res,
                const fmpz_mod_poly_t poly, ulong x, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_fit_length(res, poly->length, ctx);
    _fmpz_mod_poly_scalar_mul_ui(res->coeffs, poly->coeffs, poly->length, x, ctx);
    _fmpz_mod_poly_set_length(res, poly->length);
    _fmpz_mod_poly_normalise(res);
}
