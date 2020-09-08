/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

void
_fmpz_mod_poly_powmod_fmpz_binexp(fmpz * res, const fmpz * poly,
                                  const fmpz_t e, const fmpz * f,
                                  slong lenf, const fmpz_t p)
{
    fmpz * T, * Q;
    fmpz_t invf;
    slong lenT, lenQ;
    slong i;

    if (lenf == 2)
    {
        fmpz_powm(res, poly, e, p);
        return;
    }

    lenT = 2 * lenf - 3;
    lenQ = lenT - lenf + 1;

    T = _fmpz_vec_init(lenT + lenQ);
    Q = T + lenT;

    fmpz_init(invf);
    fmpz_invmod(invf, f + lenf - 1, p);

    _fmpz_vec_set(res, poly, lenf - 1);

    for (i = fmpz_sizeinbase(e, 2) - 2; i >= 0; i--)
    {
        _fmpz_mod_poly_sqr(T, res, lenf - 1, p);
        _fmpz_mod_poly_divrem(Q, res, T, 2 * lenf - 3, f, lenf, invf, p);

        if (fmpz_tstbit(e, i))
        {
            _fmpz_mod_poly_mul(T, res, lenf - 1, poly, lenf - 1, p);
            _fmpz_mod_poly_divrem(Q, res, T, 2 * lenf - 3, f, lenf, invf, p);
        }
    }

    fmpz_clear(invf);
    _fmpz_vec_clear(T, lenT + lenQ);
}


void
fmpz_mod_poly_powmod_fmpz_binexp(fmpz_mod_poly_t res,
                           const fmpz_mod_poly_t poly, const fmpz_t e,
                             const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
{
    fmpz * q;
    slong len = poly->length;
    slong lenf = f->length;
    slong trunc = lenf - 1;
    int qcopy = 0;

    if (lenf == 0)
    {
        flint_printf("Exception (fmpz_mod_poly_powmod_fmpz_binexp). Divide by zero\n");
        flint_abort();
    }

    if (lenf == 1)
    {
        fmpz_mod_poly_zero(res, ctx);
        return;
    }

    if (fmpz_sgn(e) < 0)
    {
        flint_printf("Exception (fmpz_mod_poly_powmod_fmpz_binexp). negative exp not implemented\n");
        flint_abort();
    }

    if (len >= lenf)
    {
        fmpz_mod_poly_t t, r;
        fmpz_mod_poly_init(t, ctx);
        fmpz_mod_poly_init(r, ctx);
        fmpz_mod_poly_divrem(t, r, poly, f, ctx);
        fmpz_mod_poly_powmod_fmpz_binexp(res, r, e, f, ctx);
        fmpz_mod_poly_clear(t, ctx);
        fmpz_mod_poly_clear(r, ctx);
        return;
    }

    if (fmpz_abs_fits_ui(e))
    {
        ulong exp = fmpz_get_ui(e);

        if (exp <= 2)
        {
            if (exp == UWORD(0))
            {
                fmpz_mod_poly_fit_length(res, 1, ctx);
                fmpz_one(res->coeffs);
                _fmpz_mod_poly_set_length(res, 1);
            }
            else if (exp == UWORD(1))
            {
                fmpz_mod_poly_set(res, poly, ctx);
            }
            else
            {
                fmpz_mod_poly_mulmod(res, poly, poly, f, ctx);
            }
            return;
        }
    }

    if (len == 0)
    {
        fmpz_mod_poly_zero(res, ctx);
        return;
    }

    if (poly->length < trunc)
    {
        q = _fmpz_vec_init(trunc);
        _fmpz_vec_set(q, poly->coeffs, len);
        _fmpz_vec_zero(q + len, trunc - len);
        qcopy = 1;
    }
    else
    {
        q = poly->coeffs;
    }

    if ((res == poly && !qcopy) || (res == f))
    {
        fmpz_mod_poly_t t;
        fmpz_mod_poly_init2(t, 2*lenf - 3, ctx);
        _fmpz_mod_poly_powmod_fmpz_binexp(t->coeffs,
                             q, e, f->coeffs, lenf, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_poly_swap(res, t, ctx);
        fmpz_mod_poly_clear(t, ctx);
    }
    else
    {
        fmpz_mod_poly_fit_length(res, 2*lenf - 3, ctx);
        _fmpz_mod_poly_powmod_fmpz_binexp(res->coeffs,
                             q, e, f->coeffs, lenf, fmpz_mod_ctx_modulus(ctx));
    }

    if (qcopy)
        _fmpz_vec_clear(q, trunc);

    _fmpz_mod_poly_set_length(res, trunc);
    _fmpz_mod_poly_normalise(res);
}
