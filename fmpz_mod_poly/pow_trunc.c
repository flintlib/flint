/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2010 William Hart
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

void
_fmpz_mod_poly_pow_trunc(fmpz * res, const fmpz * poly,
                         ulong e, slong trunc, const fmpz_t p)
{
    _fmpz_mod_poly_pow_trunc_binexp(res, poly, e, trunc, p);
}

void
fmpz_mod_poly_pow_trunc(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly,
                                ulong e, slong trunc, const fmpz_mod_ctx_t ctx)
{
    const slong len = poly->length;
    fmpz * q;
    int qcopy = 0;

    if (len < 2 || e < UWORD(3) || trunc == 0)
    {
        if (len == 0 || trunc == 0)
        {
            fmpz_mod_poly_zero(res, ctx);
        }
        else if (len == 1)
        {
            fmpz_mod_poly_fit_length(res, 1, ctx);
            fmpz_powm_ui(res->coeffs, poly->coeffs, e, fmpz_mod_ctx_modulus(ctx));
            _fmpz_mod_poly_set_length(res, 1);
            _fmpz_mod_poly_normalise(res);
        }
        else if (e == UWORD(0))
        {
            fmpz_mod_poly_set_coeff_ui(res, 0, UWORD(1), ctx);
            _fmpz_mod_poly_set_length(res, 1);
            _fmpz_mod_poly_normalise(res);
        }
        else if (e == UWORD(1))
        {
            fmpz_mod_poly_set(res, poly, ctx);
            fmpz_mod_poly_truncate(res, trunc, ctx);
        }
        else  /* e == UWORD(2) */
            fmpz_mod_poly_mullow(res, poly, poly, trunc, ctx);

        return;
    }

    if (poly->length < trunc)
    {
        q = _fmpz_vec_init(trunc);
        _fmpz_vec_set(q, poly->coeffs, poly->length);
        _fmpz_vec_zero(q + poly->length, trunc - poly->length);
        qcopy = 1;
    }
    else
    {
        q = poly->coeffs;
    }

    if (res != poly || qcopy)
    {
        fmpz_mod_poly_fit_length(res, trunc, ctx);
        _fmpz_mod_poly_pow_trunc(res->coeffs, q, e, trunc, fmpz_mod_ctx_modulus(ctx));
    }
    else
    {
        fmpz_mod_poly_t t;
        fmpz_mod_poly_init2(t, trunc, ctx);
        _fmpz_mod_poly_pow_trunc(t->coeffs, q, e, trunc, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_poly_swap(res, t, ctx);
        fmpz_mod_poly_clear(t, ctx);
    }

    if (qcopy)
        _fmpz_vec_clear(q, trunc);

    _fmpz_mod_poly_set_length(res, trunc);
    _fmpz_mod_poly_normalise(res);
}

