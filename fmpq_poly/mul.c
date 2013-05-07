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

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void _fmpq_poly_mul(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly1, const fmpz_t den1, long len1, 
                    const fmpz * poly2, const fmpz_t den2, long len2)
{
    fmpz_t gcd1;  /* GCD( poly1, den2 ) */
    fmpz_t gcd2;  /* GCD( poly2, den1 ) */
    fmpz_init(gcd1);
    fmpz_init(gcd2);
    fmpz_one(gcd1);
    fmpz_one(gcd2);
    
    if (*den2 != 1L)
    {
        _fmpz_vec_content(gcd1, poly1, len1);
        if (*gcd1 != 1L)
            fmpz_gcd(gcd1, gcd1, den2);
    }
    if (*den1 != 1L)
    {
        _fmpz_vec_content(gcd2, poly2, len2);
        if (*gcd2 != 1L)
            fmpz_gcd(gcd2, gcd2, den1);
    }
    
    /*
       TODO:  If gcd1 and gcd2 are very large compared to the degrees of 
       poly1 and poly2, we might want to create copies of the polynomials 
       and divide out the common factors *before* the multiplication.
     */
    
    _fmpz_poly_mul(rpoly, poly1, len1, poly2, len2);
    fmpz_mul(rden, den1, den2);

    if ((*gcd1 != 1L) | (*gcd2 != 1L))
    {
        fmpz_t g;
        fmpz_init(g);
        fmpz_mul(g, gcd1, gcd2);
        _fmpz_vec_scalar_divexact_fmpz(rpoly, rpoly, len1 + len2 - 1, g);
        fmpz_divexact(rden, rden, g);
        fmpz_clear(g);
    }
    
    fmpz_clear(gcd1);
    fmpz_clear(gcd2);
}

void fmpq_poly_mul(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
{
    long len;
    
    if (poly1->length == 0 || poly2->length == 0)
    {
        fmpq_poly_zero(res);
        return;
    }

    len = poly1->length + poly2->length - 1;

    if (res == poly2 || res == poly1)
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, len);
        fmpq_poly_mul(t, poly1, poly2);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
        return;
    }
    
    fmpq_poly_fit_length(res, len);

    if (poly1->length >= poly2->length)
        _fmpq_poly_mul(res->coeffs, res->den, 
                       poly1->coeffs, poly1->den, poly1->length,
                       poly2->coeffs, poly2->den, poly2->length);
    else
        _fmpq_poly_mul(res->coeffs, res->den, 
                       poly2->coeffs, poly2->den, poly2->length,
                       poly1->coeffs, poly1->den, poly1->length);

    _fmpq_poly_set_length(res, len);
}

