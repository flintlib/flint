/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

void
_fmpz_mod_poly_powmod_fmpz_binexp(fmpz * res, const fmpz * poly,
                                  const fmpz_t e, const fmpz * f,
                                  len_t lenf, const fmpz_t p)
{
    fmpz * T, * Q;
    fmpz_t invf;
    len_t lenT, lenQ;
    len_t i;

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
                           const fmpz_mod_poly_t f)
{
    fmpz * q;
    len_t len = poly->length;
    len_t lenf = f->length;
    len_t trunc = lenf - 1;
    int qcopy = 0;

    if (lenf == 0)
    {
        printf("Exception: fmpz_mod_poly_powmod: divide by zero\n");
        abort();
    }

    if (fmpz_sgn(e) < 0)
    {
        printf("Exception: fmpz_mod_poly_powmod: negative exp not implemented\n");
        abort();
    }

    if (len >= lenf)
    {
        fmpz_mod_poly_t t, r;
        fmpz_mod_poly_init(t, &res->p);
        fmpz_mod_poly_init(r, &res->p);
        fmpz_mod_poly_divrem(t, r, poly, f);
        fmpz_mod_poly_powmod_fmpz_binexp(res, r, e, f);
        fmpz_mod_poly_clear(t);
        fmpz_mod_poly_clear(r);
        return;
    }

    if (fmpz_abs_fits_ui(e))
    {
        ulong exp = fmpz_get_ui(e);

        if (exp <= 2)
        {
            if (exp == 0UL)
            {
                fmpz_mod_poly_fit_length(res, 1);
                fmpz_one(res->coeffs);
                _fmpz_mod_poly_set_length(res, 1);
            }
            else if (exp == 1UL)
            {
                fmpz_mod_poly_set(res, poly);
            }
            else
                fmpz_mod_poly_mulmod(res, poly, poly, f);
            return;
        }
    }

    if (lenf == 1 || len == 0)
    {
        fmpz_mod_poly_zero(res);
        return;
    }

    if (poly->length < trunc)
    {
        q = _fmpz_vec_init(trunc);
        _fmpz_vec_set(q, poly->coeffs, len);
        _fmpz_vec_zero(q + len, trunc - len);
        qcopy = 1;
    } else
        q = poly->coeffs;

    if ((res == poly && !qcopy) || (res == f))
    {
        fmpz_mod_poly_t t;
        fmpz_mod_poly_init2(t, &poly->p, 2 * lenf - 3);
        _fmpz_mod_poly_powmod_fmpz_binexp(t->coeffs,
            q, e, f->coeffs, lenf, &poly->p);
        fmpz_mod_poly_swap(res, t);
        fmpz_mod_poly_clear(t);
    }
    else
    {
        fmpz_mod_poly_fit_length(res, 2 * lenf - 3);
        _fmpz_mod_poly_powmod_fmpz_binexp(res->coeffs,
            q, e, f->coeffs, lenf, &poly->p);
    }

    if (qcopy)
        _fmpz_vec_clear(q, trunc);

    _fmpz_mod_poly_set_length(res, trunc);
    _fmpz_mod_poly_normalise(res);
}
