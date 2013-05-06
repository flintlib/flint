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
_fmpq_poly_compose_series(fmpz * res, fmpz_t den, const fmpz * poly1,
        const fmpz_t den1, long len1, const fmpz * poly2,
        const fmpz_t den2, long len2, long n)
{
    if (len1 <= 20)
        _fmpq_poly_compose_series_horner(res, den, poly1, den1, len1,
                poly2, den2, len2, n);
    else
        _fmpq_poly_compose_series_brent_kung(res, den, poly1, den1, len1,
                poly2, den2, len2, n);
}

void
fmpq_poly_compose_series(fmpq_poly_t res, 
                    const fmpq_poly_t poly1, const fmpq_poly_t poly2, long n)
{
    long len1 = poly1->length;
    long len2 = poly2->length;
    long lenr;

    if (len2 != 0 && !fmpz_is_zero(poly2->coeffs))
    {
        printf("Exception (fmpq_poly_compose_series). Inner polynomial \n"
               "must have zero constant term,\n");
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
        _fmpq_poly_compose_series(res->coeffs, res->den, 
                           poly1->coeffs, poly1->den, len1, 
                           poly2->coeffs, poly2->den, len2, lenr);
        _fmpq_poly_set_length(res, lenr);
        _fmpq_poly_normalise(res);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, lenr);
        _fmpq_poly_compose_series(t->coeffs, t->den, 
                           poly1->coeffs, poly1->den, len1,
                           poly2->coeffs, poly2->den, len2, lenr);
        _fmpq_poly_set_length(t, lenr);
        _fmpq_poly_normalise(t);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }
}
