/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#undef ulong
#define ulong ulongxx/* interferes with system includes */

#include <stdlib.h>

#undef ulong

#include <gmp.h>

#define ulong mp_limb_t

#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"
#include "long_extras.h"

void
_fmpz_mod_poly_powmod_x_fmpz_preinv(fmpz * res, const fmpz_t e, const fmpz * f,
                                    slong lenf, const fmpz* finv, slong lenfinv,
                                    const fmpz_t p)
{
    fmpz * T, * Q;
    slong lenT, lenQ;
    slong i, window, l, c;

    lenT = 2 * lenf - 3;
    lenQ = lenT - lenf + 1;

    T = _fmpz_vec_init(lenT + lenQ);
    Q = T + lenT;

    fmpz_one(res);
    _fmpz_vec_zero(res + 1, lenf - 2);
    l = z_sizeinbase(lenf - 1, 2) - 2;
    window = (WORD(1) << l);
    c = l;
    i = fmpz_sizeinbase(e, 2) - 2;
    if (i <= l)
    {
      window = (WORD(1) << i);
      c = i;
      l = i;
    }

    if (c == 0)
    {
        _fmpz_mod_poly_shift_left(T, res, lenf - 1, window);
        _fmpz_mod_poly_divrem_newton_n_preinv(Q, res, T, lenf - 1 + window, f,
                                              lenf, finv, lenfinv, p);
        c = l + 1;
        window = WORD(0);
    }

    for (; i >= 0; i--)
    {
        _fmpz_mod_poly_sqr(T, res, lenf - 1, p);
        _fmpz_mod_poly_divrem_newton_n_preinv(Q, res, T, 2 * lenf - 3, f, lenf,
                                              finv, lenfinv, p);

        c--;
        if (fmpz_tstbit(e, i))
        {
            if (window == WORD(0) && i <= l - 1)
                c = i;
            if ( c >= 0)
              window = window | (WORD(1) << c);
        }
        else if (window == WORD(0))
            c = l + 1;
        if (c == 0)
        {
            _fmpz_mod_poly_shift_left(T, res, lenf - 1, window);
            
            _fmpz_mod_poly_divrem_newton_n_preinv(Q, res, T, lenf - 1 + window,
                                                  f, lenf, finv, lenfinv, p);
            c = l + 1;
            window = WORD(0);
        }
    }

    _fmpz_vec_clear(T, lenT + lenQ);
}


void
fmpz_mod_poly_powmod_x_fmpz_preinv(fmpz_mod_poly_t res, const fmpz_t e,
                           const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv,
                                                      const fmpz_mod_ctx_t ctx)
{
    slong lenf = f->length;
    slong trunc = lenf - 1;
    fmpz_mod_poly_t tmp;

    if (lenf == 0)
    {
        flint_printf("Exception (fmpz_mod_poly_powmod_x_fmpz_preinv)."
                     "Divide by zero\n");
        flint_abort();
    }

    if (fmpz_sgn(e) < 0)
    {
        flint_printf("Exception (fmpz_mod_poly_powmod_x_fmpz_preinv)."
                     "Negative exp not implemented\n");
        flint_abort();
    }

    if (lenf == 1)
    {
        fmpz_mod_poly_zero(res, ctx);
        return;
    }

    if (lenf == 2)
    {
        fmpz_mod_poly_t r, poly;
        fmpz_mod_poly_init(tmp, ctx);
        fmpz_mod_poly_init(r, ctx);
        fmpz_mod_poly_init2(poly, 2, ctx);
        fmpz_mod_poly_set_coeff_ui(poly, 1, 1, ctx);
        fmpz_mod_poly_divrem(tmp, r, poly, f, ctx);
        fmpz_mod_poly_powmod_fmpz_binexp_preinv(res, r, e, f, finv, ctx);
        fmpz_mod_poly_clear(tmp, ctx);
        fmpz_mod_poly_clear(r, ctx);
        fmpz_mod_poly_clear(poly, ctx);
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
                fmpz_mod_poly_t r;
                fmpz_mod_poly_init2(r, 2, ctx);
                fmpz_mod_poly_set_coeff_ui(r, 1, 1, ctx);
                fmpz_mod_poly_init(tmp, ctx);
                fmpz_mod_poly_divrem(tmp, res, r, f, ctx);
                fmpz_mod_poly_clear(tmp, ctx);
                fmpz_mod_poly_clear(r, ctx);
            }
            else
            {
                fmpz_mod_poly_init2(tmp, 3, ctx);
                fmpz_mod_poly_set_coeff_ui(tmp, 1, 1, ctx);
                fmpz_mod_poly_mulmod(res, tmp, tmp, f, ctx);
                fmpz_mod_poly_clear(tmp, ctx);
            }
            return;
        }
    }

    if ((res == f) || (res == finv))
    {
        fmpz_mod_poly_init2(tmp, trunc, ctx);
        _fmpz_mod_poly_powmod_x_fmpz_preinv(tmp->coeffs, e, f->coeffs, lenf,
                        finv->coeffs, finv->length, fmpz_mod_ctx_modulus(ctx));
        fmpz_mod_poly_swap(res, tmp, ctx);
        fmpz_mod_poly_clear(tmp, ctx);
    }
    else
    {
        fmpz_mod_poly_fit_length(res, trunc, ctx);
        _fmpz_mod_poly_powmod_x_fmpz_preinv(res->coeffs, e, f->coeffs, lenf,
                       finv->coeffs, finv->length, fmpz_mod_ctx_modulus(ctx));
    }

    _fmpz_mod_poly_set_length(res, trunc);
    _fmpz_mod_poly_normalise(res);
}
