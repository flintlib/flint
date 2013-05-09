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

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void
_fmpq_poly_mullow(fmpz * rpoly, fmpz_t rden, 
                  const fmpz * poly1, const fmpz_t den1, len_t len1, 
                  const fmpz * poly2, const fmpz_t den2, len_t len2, len_t n)
{
    _fmpz_poly_mullow(rpoly, poly1, len1, poly2, len2, n);
    fmpz_mul(rden, den1, den2);
}

void
fmpq_poly_mullow(fmpq_poly_t res,
                 const fmpq_poly_t poly1, const fmpq_poly_t poly2, len_t n)
{
    const len_t len1 = poly1->length;
    const len_t len2 = poly2->length;
    len_t lenr;

    if (len1 == 0 || len2 == 0 || n == 0)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        fmpq_poly_mullow(t, poly1, poly2, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
        return;
    }

    lenr = len1 + len2 - 1;
    if (n > lenr)
        n = lenr;

    fmpq_poly_fit_length(res, n);
    if (len1 >= len2)
        _fmpq_poly_mullow(res->coeffs, res->den, 
                          poly1->coeffs, poly1->den, len1, 
                          poly2->coeffs, poly2->den, len2, n);
    else
        _fmpq_poly_mullow(res->coeffs, res->den, 
                          poly2->coeffs, poly2->den, len2, 
                          poly1->coeffs, poly1->den, len1, n);
    _fmpq_poly_set_length(res, n);
    fmpq_poly_canonicalise(res);
}
