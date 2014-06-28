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
    Copyright (C) 2013 Martin Lee
    Copyright (C) 2014 William Hart

******************************************************************************/

#undef ulong
#define ulong ulongxx/* interferes with system includes */

#include <stdlib.h>

#undef ulong

#include <gmp.h>

#define ulong mp_limb_t

#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

void
_fmpz_mod_poly_powmod_fmpz_binexp_powers_precomp(fmpz * res, const fmpz * poly,
                                  const fmpz_t e, const fmpz * f,
                                  slong lenf, fmpz ** const powers,
                                  const fmpz_t p)
{
    fmpz * T;
    slong lenT;
    slong i;

    if (lenf == 2)
    {
        fmpz_powm(res, poly, e, p);
        return;
    }

    lenT = 2 * lenf - 3;
    
    T = _fmpz_vec_init(lenT);
    
    _fmpz_vec_set(res, poly, lenf - 1);

    for (i = fmpz_sizeinbase(e, 2) - 2; i >= 0; i--)
    {
        _fmpz_mod_poly_sqr(T, res, lenf - 1, p);
        _fmpz_mod_poly_rem_powers_precomp(T, 2 * lenf - 3, f, lenf,
                                              powers, p);
        _fmpz_vec_set(res, T, lenf - 1);
        

        if (fmpz_tstbit(e, i))
        {
            _fmpz_mod_poly_mul(T, res, lenf - 1, poly, lenf - 1, p);
            _fmpz_mod_poly_rem_powers_precomp(T, 2 * lenf - 3, f,
                                                  lenf, powers, p);
            _fmpz_vec_set(res, T, lenf - 1);
        }
    }

    _fmpz_vec_clear(T, lenT);
}


void
fmpz_mod_poly_powmod_fmpz_binexp_powers_precomp(fmpz_mod_poly_t res,
                           const fmpz_mod_poly_t poly, const fmpz_t e,
              const fmpz_mod_poly_t f, const fmpz_mod_poly_powers_precomp_t finv)
{
    fmpz * q;
    slong len = poly->length;
    slong lenf = f->length;
    slong trunc = lenf - 1;
    int qcopy = 0;

    if (lenf == 0)
    {
        flint_printf("Exception (fmpz_mod_poly_powmod_fmpz_binexp_powers_precomp)."
                     "Divide by zero.\n");
        abort();
    }

    if (fmpz_sgn(e) < 0)
    {
        flint_printf("Exception (fmpz_mod_poly_powmod_fmpz_binexp_powers_precomp)."
                     "Negative exp not implemented\n");
        abort();
    }

    if (len >= lenf)
    {
        fmpz_mod_poly_t t, r;
        fmpz_mod_poly_init(t, &res->p);
        fmpz_mod_poly_init(r, &res->p);
        fmpz_mod_poly_divrem(t, r, poly, f);
        fmpz_mod_poly_powmod_fmpz_binexp_powers_precomp(res, r, e, f, finv);
        fmpz_mod_poly_clear(t);
        fmpz_mod_poly_clear(r);
        return;
    }

    if (fmpz_abs_fits_ui(e))
    {
        ulong exp = fmpz_get_ui(e);

        if (exp <= 2)
        {
            if (exp == UWORD (0))
            {
                fmpz_mod_poly_fit_length(res, 1);
                fmpz_one(res->coeffs);
                _fmpz_mod_poly_set_length(res, 1);
            }
            else if (exp == UWORD (1))
            {
                fmpz_mod_poly_set(res, poly);
            }
            else
            {
                fmpz_mod_poly_t temp;
                fmpz_mod_poly_init(temp, &poly->p);
                fmpz_mod_poly_mul(temp, poly, poly);
                fmpz_mod_poly_rem_powers_precomp(res, temp, f, finv);
                fmpz_mod_poly_clear(temp);
            }
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
        _fmpz_mod_poly_powmod_fmpz_binexp_powers_precomp(t->coeffs,
            q, e, f->coeffs, lenf, finv->powers, &poly->p);
        fmpz_mod_poly_swap(res, t);
        fmpz_mod_poly_clear(t);
    }
    else
    {
        fmpz_mod_poly_fit_length(res, 2 * lenf - 3);
        _fmpz_mod_poly_powmod_fmpz_binexp_powers_precomp(res->coeffs,
            q, e, f->coeffs, lenf, finv->powers, &poly->p);
    }

    if (qcopy)
        _fmpz_vec_clear(q, trunc);

    _fmpz_mod_poly_set_length(res, trunc);
    _fmpz_mod_poly_normalise(res);
}
