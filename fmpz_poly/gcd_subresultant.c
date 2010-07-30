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

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void _fmpz_poly_gcd_subresultant(fmpz * res, const fmpz * poly1, 
                                 long len1, const fmpz * poly2, long len2)
{
    fmpz_t a, b, d;
    fmpz * A, * B, * Q, * R, * Z;
    long lenA, lenB, lenQ, lenR, lenZ;
    long s;
    
    fmpz * g;
    fmpz_t h, one, r, temp;
    long olddelta;
    
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(d);
    _fmpz_poly_content(a, poly1, len1);
    _fmpz_poly_content(b, poly2, len2);
    fmpz_gcd(d, a, b);
    
    if (len2 == 1L)
    {
        fmpz_set(res, d);
        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(d);
        return;
    }
    
    lenA = len1;
    lenB = len2;
    lenQ = lenA + lenB - 1L;
    lenR = lenA;
    lenZ = lenA + lenB + lenQ + lenR;
    Z = _fmpz_vec_init(lenZ);
    A = Z;
    B = A + lenA;
    Q = B + lenB;
    R = Q + lenQ;
    _fmpz_vec_scalar_divexact_fmpz(A, poly1, len1, a);
    _fmpz_vec_scalar_divexact_fmpz(B, poly2, len2, b);
    s = 0L;
    
    fmpz_init(h);
    fmpz_init(one);
    fmpz_init(r);
    fmpz_init(temp);
    fmpz_set_ui(h, 1UL);
    fmpz_set_ui(one, 1UL);
    g = one;
    olddelta = 1L;
    
    while (1)
    {
        long delta = lenA - lenB;
        _fmpz_poly_pseudo_divrem(Q, R, &s, A, lenA, B, lenB);
        lenQ = lenA + (lenB - 1L);
        lenR = lenB - 1L;
        for ( ; lenR != 0L && fmpz_is_zero(R + (lenR - 1L)); lenR--) ;
        if (lenR <= 1L)
            break;
        
        {   /* Swap A and B */
            fmpz * T;
            long t;
            T = A, A = B, B = T, t = lenA, lenA = lenB, lenB = t;
        }
        
        if (olddelta == 1L)
            fmpz_pow_ui(r, g, delta + 1L);
        else
        {
            fmpz_pow_ui(r, h, delta);
            fmpz_mul(r, r, g);
        }
        
        g = A + (lenA - 1L);
        fmpz_pow_ui(temp, g, delta + 1L - s);
        _fmpz_vec_scalar_mul_fmpz(R, R, lenR, temp);
        
        _fmpz_vec_scalar_divexact_fmpz(B, R, lenR, r);
        lenB = lenR;
        
        olddelta = delta;
        if (delta == 0L)
            fmpz_set_ui(h, 1UL);
        else if (delta == 1L)
            fmpz_set(h, g);
        else
        {
            fmpz_t num, den;
            fmpz_init(num), fmpz_init(den);
            fmpz_pow_ui(num, g, delta);
            fmpz_pow_ui(den, h, delta - 1UL);
            fmpz_divexact(h, num, den);
            fmpz_clear(num), fmpz_clear(den);
        }
    }
    if (lenR == 1L)
    {
        fmpz_set_ui(B, 1UL);
        lenB = 1L;
    }
    
    /* {res, len2} = +- (d / b) {B, lenB} */
    _fmpz_poly_content(b, B, lenB);
    _fmpz_vec_scalar_divexact_fmpz(res, B, lenB, b);
    if (fmpz_sgn(res + (lenB - 1L)) < 0)
        fmpz_neg(d, d);
    _fmpz_vec_scalar_mul_fmpz(res, res, lenB, d);
    for (s = lenB; s < len2; s++)
        fmpz_zero(res + s);
    
    _fmpz_vec_clear(Z, lenZ);
    fmpz_clear(a);
    fmpz_clear(b);
    fmpz_clear(d);
    fmpz_clear(h);
    fmpz_clear(one);
    fmpz_clear(r);
    fmpz_clear(temp);
}   

void fmpz_poly_gcd_subresultant(fmpz_poly_t res, 
                    const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    long len1 = poly1->length;
    long len2 = poly2->length;
    
    if (len1 == 0)
    {
        if (len2 == 0)
            fmpz_poly_zero(res);
        else
        {
            if (fmpz_sgn(poly2->coeffs + (len2 - 1)) > 0)
                fmpz_poly_set(res, poly2);
            else
                fmpz_poly_neg(res, poly2);
        }
        return;
    }
    else
    {
        if (len2 == 0)
        {
            if (fmpz_sgn(poly1->coeffs + (len1 - 1)) > 0)
                fmpz_poly_set(res, poly1);
            else
                fmpz_poly_neg(res, poly1);
            return;
        }
    }
    
    long rlen = FLINT_MIN(len1, len2);
    fmpz_poly_fit_length(res, rlen);
    if (len1 >= len2)
        _fmpz_poly_gcd_subresultant(res->coeffs, poly1->coeffs, len1, 
            poly2->coeffs, len2);
    else 
        _fmpz_poly_gcd_subresultant(res->coeffs, poly2->coeffs, len2, 
            poly1->coeffs, len1);
    _fmpz_poly_set_length(res, rlen);
    _fmpz_poly_normalise(res);
}
