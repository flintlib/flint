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
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void
_fmpq_poly_compose_series_horner(fmpz * res, fmpz_t den, const fmpz * poly1,
        const fmpz_t den1, long len1, const fmpz * poly2,
        const fmpz_t den2, long len2, long n)
{
    if (fmpz_is_one(den2))
    {
        _fmpz_poly_compose_series(res, poly1, len1, poly2, len2, n);
        fmpz_set(den, den1);
        _fmpq_poly_canonicalise(res, den, n);
    }
    else if (n == 1)
    {
        fmpz_set(res, poly1);
        fmpz_set(den, den1);
        _fmpq_poly_canonicalise(res, den, 1);
    }
    else
    {
        long i = len1 - 1;
        long lenr;
        fmpz_t tden;
        fmpz * t = _fmpz_vec_init(n);
        fmpz_init(tden);

        _fmpz_vec_zero(res, n);

        lenr = len2;
        _fmpq_poly_scalar_mul_fmpz(res, den, poly2, den2, len2, poly1 + i);
        _fmpq_poly_scalar_div_fmpz(res, den, res, den, len2, den1);
        i--;

        _fmpq_poly_add(res, den, res, den, len2, poly1 + i, den1, 1);
        _fmpq_poly_canonicalise(res, den, lenr);

        while (i > 0)
        {
            i--;
            if (lenr + len2 - 1 < n)
            {
                _fmpq_poly_mul(t, tden, res, den, lenr, poly2, den2, len2);
                lenr = lenr + len2 - 1;
            }
            else
            {
                _fmpq_poly_mullow(t, tden, res, den, lenr,
                                            poly2, den2, len2, n);
                lenr = n;
            }
            _fmpq_poly_canonicalise(t, tden, lenr);
            _fmpq_poly_add(res, den, t, tden, lenr, poly1 + i, den1, 1);
        }

        _fmpq_poly_canonicalise(res, den, n);

        _fmpz_vec_clear(t, n);
        fmpz_clear(tden);
    }
}

void
fmpq_poly_compose_series_horner(fmpq_poly_t res, 
                    const fmpq_poly_t poly1, const fmpq_poly_t poly2, long n)
{
    long len1 = poly1->length;
    long len2 = poly2->length;
    long lenr;

    if (len2 != 0 && !fmpz_is_zero(poly2->coeffs))
    {
        printf("Exception (fmpq_poly_compose_series_horner). Inner polynomial \n"
               "must have zero constant term.\n");
        abort();
    }

    if (len1 == 0 || n == 0)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (len2 == 0 || len1 == 1)
    {
        fmpq_poly_fit_length(res, 1);
        fmpz_set(res->coeffs, poly1->coeffs);
        fmpz_set(res->den, poly1->den);
        {
            fmpz_t d;
            fmpz_init(d);
            fmpz_gcd(d, res->coeffs, res->den);
            if (!fmpz_is_one(d))
            {
                fmpz_divexact(res->coeffs, res->coeffs, d);
                fmpz_divexact(res->den, res->den, d);
            }
            fmpz_clear(d);
        }
        _fmpq_poly_set_length(res, 1);
        _fmpq_poly_normalise(res);
        return;
    }

    lenr = FLINT_MIN((len1 - 1) * (len2 - 1) + 1, n);
    len1 = FLINT_MIN(len1, lenr);
    len2 = FLINT_MIN(len2, lenr);

    if ((res != poly1) && (res != poly2))
    {
        fmpq_poly_fit_length(res, lenr);
        _fmpq_poly_compose_series_horner(res->coeffs, res->den, 
                           poly1->coeffs, poly1->den, len1, 
                           poly2->coeffs, poly2->den, len2, lenr);
        _fmpq_poly_set_length(res, lenr);
        _fmpq_poly_normalise(res);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, lenr);
        _fmpq_poly_compose_series_horner(t->coeffs, t->den, 
                           poly1->coeffs, poly1->den, len1,
                           poly2->coeffs, poly2->den, len2, lenr);
        _fmpq_poly_set_length(t, lenr);
        _fmpq_poly_normalise(t);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }
}
