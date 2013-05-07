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

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

void
_fmpz_mod_poly_compose_mod_horner(fmpz * res, const fmpz * f, long lenf, const fmpz * g,
                                              const fmpz * h, long lenh, const fmpz_t p)
{
    long i, len;
    fmpz * t;

    if (lenh == 1)
        return;

    if (lenf == 1)
    {
        fmpz_set(res, f);
        return;
    }

    if (lenh == 2)
    {
        _fmpz_mod_poly_evaluate_fmpz(res, f, lenf, g, p);
        return;
    }

    len = lenh - 1;
    i = lenf - 1;
    t = _fmpz_vec_init(2 * lenh - 3);

    _fmpz_mod_poly_scalar_mul_fmpz(res, g, len, f + i, p);
    i--;
    if (i >= 0)
    {
        fmpz_add(res, res, f + i);
        fmpz_mod(res, res, p);
    }

    while (i > 0)
    {
        i--;
        _fmpz_mod_poly_mulmod(t, res, len, g, len, h, lenh, p);
        _fmpz_mod_poly_add(res, t, len, f + i, 1, p);
    }

    _fmpz_vec_clear(t, 2 * lenh - 3);
}

void
fmpz_mod_poly_compose_mod_horner(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1,
                         const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t poly3)
{
    fmpz_t inv3;
    long len1 = poly1->length;
    long len2 = poly2->length;
    long len3 = poly3->length;
    long len = len3 - 1;
    long vec_len = FLINT_MAX(len3 - 1, len2);

    fmpz * ptr2;

    if (len3 == 0)
    {
        printf("Exception: division by zero in fmpz_mod_poly_compose_mod_horner\n");
        abort();
    }

    if (len1 == 0 || len3 == 1)
    {
        fmpz_mod_poly_zero(res);
        return;
    }

    if (len1 == 1)
    {
        fmpz_mod_poly_set(res, poly1);
        return;
    }

    if (res == poly3 || res == poly1)
    {
        fmpz_mod_poly_t tmp;
        fmpz_mod_poly_init(tmp, &res->p);
        fmpz_mod_poly_compose_mod_horner(tmp, poly1, poly2, poly3);
        fmpz_mod_poly_swap(tmp, res);
        fmpz_mod_poly_clear(tmp);
        return;
    }

    ptr2 = _fmpz_vec_init(vec_len);

    if (len2 <= len3 - 1)
    {
        _fmpz_vec_set(ptr2, poly2->coeffs, len2);
        _fmpz_vec_zero(ptr2 + len2, vec_len - len2);
    }
    else
    {
        fmpz_init(inv3);
        fmpz_invmod(inv3, poly3->coeffs + len, &res->p);
        _fmpz_mod_poly_rem(ptr2, poly2->coeffs, len2,
                           poly3->coeffs, len3, inv3, &res->p);
        fmpz_clear(inv3);
    }

    fmpz_mod_poly_fit_length(res, len);
    _fmpz_mod_poly_compose_mod_horner(res->coeffs,
        poly1->coeffs, len1, ptr2, poly3->coeffs, len3, &res->p);
    _fmpz_mod_poly_set_length(res, len);
    _fmpz_mod_poly_normalise(res);

    _fmpz_vec_clear(ptr2, vec_len);
}
