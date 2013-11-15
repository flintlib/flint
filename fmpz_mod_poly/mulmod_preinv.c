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

    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Martin Lee

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

void _fmpz_mod_poly_mulmod_preinv(fmpz * res, const fmpz * poly1, slong len1,
                    const fmpz * poly2, slong len2, const fmpz * f, slong lenf,
                    const fmpz* finv, slong lenfinv, const fmpz_t p)
{
    fmpz * T, * Q;
    slong lenT, lenQ;

    lenT = len1 + len2 - 1;
    lenQ = lenT - lenf + 1;

    T = _fmpz_vec_init(lenT + lenQ);
    Q = T + lenT;

    if (len1 >= len2)
        _fmpz_mod_poly_mul(T, poly1, len1, poly2, len2, p);
    else
        _fmpz_mod_poly_mul(T, poly2, len2, poly1, len1, p);

    _fmpz_mod_poly_divrem_newton_n_preinv(Q, res, T, lenT, f, lenf,
                                          finv, lenfinv, p);

    _fmpz_vec_clear(T, lenT + lenQ);
}

void
fmpz_mod_poly_mulmod_preinv(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                         const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t f,
                         const fmpz_mod_poly_t finv)
{
    slong len1, len2, lenf;
    fmpz * fcoeffs;

    lenf = f->length;
    len1 = poly1->length;
    len2 = poly2->length;

    if (lenf == 0)
    {
        flint_printf("Exception (fmpz_mod_poly_mulmod_preinv). Divide by zero\n");
        abort();
    }

    if (lenf <= len1 || lenf <= len2)
    {
        flint_printf("Exception (fmpz_mod_poly_mulmod_preinv). Input larger than modulus.\n");
        abort();
    }

    if (lenf == 1 || len1 == 0 || len2 == 0)
    {
        fmpz_mod_poly_zero(res);
        return;
    }

    if (len1 + len2 - lenf > 0)
    {
        if (f == res)
        {
            fcoeffs = _fmpz_vec_init(lenf);
            _fmpz_vec_set(fcoeffs, f->coeffs, lenf);
        }
        else
            fcoeffs = f->coeffs;

        fmpz_mod_poly_fit_length(res, len1 + len2 - 1);
        _fmpz_mod_poly_mulmod_preinv(res->coeffs, poly1->coeffs, len1,
                              poly2->coeffs, len2, fcoeffs, lenf,
                              finv->coeffs, finv->length, &res->p);
        if (f == res)
            _fmpz_vec_clear(fcoeffs, lenf);

        _fmpz_mod_poly_set_length(res, lenf - 1);
        _fmpz_mod_poly_normalise(res);
    }
    else
    {
        fmpz_mod_poly_mul(res, poly1, poly2);
    }
}
