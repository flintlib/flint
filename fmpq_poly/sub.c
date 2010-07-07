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
    Copyright (C) 2008, 2009 William Hart

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void _fmpq_poly_sub_in_place(fmpz * rpoly, fmpz_t rden, ulong rlen, 
                          const fmpz * poly, const fmpz_t den, const ulong len)
{
    ulong i;
    
    if (*rden == 1)
    {
        if (*den == 1)
            _fmpz_vec_sub(rpoly, rpoly, poly, len);
        else
        {
            _fmpz_vec_scalar_mul_fmpz(rpoly, rpoly, rlen, den);
            _fmpz_vec_sub(rpoly, rpoly, poly, len);
            fmpz_set(rden, den);
        }
    }
    else
    {
        if (*den == 1)
            _fmpz_vec_scalar_submul_fmpz(rpoly, poly, len, rden);
        else
        {
            fmpz_t d;
            fmpz_init(d);
            fmpz_gcd(d, rden, den);
            if (*d == 1)
            {
                /* Set rpoly = den rpoly - rden poly, rden = rden den. */
                _fmpz_vec_scalar_mul_fmpz(rpoly, rpoly, rlen, den);
                _fmpz_vec_scalar_submul_fmpz(rpoly, poly, len, rden);
                fmpz_mul(rden, rden, den);
            }
            else
            {
                fmpz_t _rden, _den;
                fmpz_init(_rden);
                fmpz_init(_den);
                fmpz_divexact(_rden, rden, d);
                fmpz_divexact(_den, den, d);
                
                /* Set rpoly = den rpoly - rden poly. */
                _fmpz_vec_scalar_mul_fmpz(rpoly, rpoly, rlen, _den);
                _fmpz_vec_scalar_submul_fmpz(rpoly, poly, len, _rden);
                
                if (_fmpz_vec_is_zero(rpoly, rlen))
                    fmpz_set_si(rden, 1);
                else
                {
                    fmpz_mul(rden, rden, _den);
                    fmpz_t content;
                    fmpz_init(content);
                    _fmpz_vec_content(content, rpoly, rlen);
                    fmpz_gcd(d, content, d);
                    if (*d != 1)
                    {
                        _fmpz_vec_scalar_divexact(rpoly, rpoly, rlen, d);
                        fmpz_divexact(rden, rden, d);
                    }
                    fmpz_clear(content);
                }
                fmpz_clear(_rden);
                fmpz_clear(_den);
            }
            fmpz_clear(d);
        }
    }
}

void _fmpq_poly_sub(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly1, const fmpz_t den1, const ulong len1, 
                    const fmpz * poly2, const fmpz_t den2, const ulong len2)
{
    ulong max = FLINT_MAX(len1, len2);
    ulong min = FLINT_MIN(len1, len2);
    ulong i;
    
    if (*den1 == 1)
    {
        if (*den2 == 1)
        {
            for (i = 0; i < min; i++)
                fmpz_sub(rpoly + i, poly1 + i, poly2 + i);
            for (i = min; i < len1; i++)
                fmpz_set(rpoly + i, poly1 + i);
            for (i = min; i < len2; i++)
                fmpz_neg(rpoly + i, poly2 + i);
            fmpz_set_si(rden, 1);
        }
        else
        {
            for (i = 0; i < len2; i++)
                fmpz_neg(rpoly + i, poly2 + i);
            for ( ; i < max; i++)
                _fmpz_demote(rpoly + i);
            _fmpz_vec_scalar_addmul_fmpz(rpoly, poly1, len1, den2);
            fmpz_set(rden, den2);
        }
    }
    else
    {
        if (*den2 == 1)
        {
            for (i = 0; i < len1; i++)
                fmpz_set(rpoly + i, poly1 + i);
            for ( ; i < max; i++)
                _fmpz_demote(rpoly + i);
            _fmpz_vec_scalar_submul_fmpz(rpoly, poly2, len2, den1);
            fmpz_set(rden, den1);
        }
        else
        {
            fmpz_t d;
            fmpz_init(d);
            fmpz_gcd(d, den1, den2);
            if (*d == 1)
            {
                /* First set rpoly = den2 poly1 - den1 poly2 up to min,      */
                /* then deal with the remaining entries in either poly1 or   */
                /* poly2.                                                    */
                _fmpz_vec_scalar_mul_fmpz(rpoly, poly1, min, den2);
                _fmpz_vec_scalar_submul_fmpz(rpoly, poly2, min, den1);
                _fmpz_vec_scalar_mul_fmpz(rpoly + min, poly1 + min, len1 - min, den2);
                _fmpz_vec_scalar_mul_fmpz(rpoly + min, poly2 + min, len2 - min, den1);
                _fmpz_vec_neg(rpoly + min, rpoly + min, len2 - min);
                fmpz_mul(rden, den1, den2);
            }
            else
            {
                fmpz_t den11, den22;
                fmpz_init(den11);
                fmpz_init(den22);
                fmpz_divexact(den11, den1, d);
                fmpz_divexact(den22, den2, d);
                
                /* First set rpoly = den22 poly1 - den11 poly2 up to min,    */
                /* then deal with the remaining entries in either poly1 or   */
                /* poly2.                                                    */
                _fmpz_vec_scalar_mul_fmpz(rpoly, poly1, min, den22);
                _fmpz_vec_scalar_submul_fmpz(rpoly, poly2, min, den11);
                _fmpz_vec_scalar_mul_fmpz(rpoly + min, poly1 + min, len1 - min, den22);
                _fmpz_vec_scalar_mul_fmpz(rpoly + min, poly2 + min, len2 - min, den11);
                _fmpz_vec_neg(rpoly + min, rpoly + min, len2 - min);
                
                if (_fmpz_vec_is_zero(rpoly, max))
                    fmpz_set_si(rden, 1);
                else
                {
                    fmpz_mul(rden, den1, den22);
                    fmpz_t content;
                    fmpz_init(content);
                    _fmpz_vec_content(content, rpoly, max);
                    fmpz_gcd(d, content, d);
                    if (*d != 1)
                    {
                        _fmpz_vec_scalar_divexact(rpoly, rpoly, max, d);
                        fmpz_divexact(rden, rden, d);
                    }
                    fmpz_clear(content);
                }
                fmpz_clear(den11);
                fmpz_clear(den22);
            }
            fmpz_clear(d);
        }
    }
}

void fmpq_poly_sub_naive(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
{
    ulong i, max, min;
    fmpz_t x, y;
    
    max = FLINT_MAX(poly1->length, poly2->length);
    min = FLINT_MIN(poly1->length, poly2->length);
    
    fmpq_poly_fit_length(res, max);
    
    fmpz_init(x);
    fmpz_init(y);
    for (i = 0; i < min; i++)
    {
        fmpz_mul(x, poly2->den, poly1->coeffs + i);
        fmpz_mul(y, poly1->den, poly2->coeffs + i);
        fmpz_sub(res->coeffs + i, x, y);
    }
    for (i = min; i < poly1->length; i++)
        fmpz_mul(res->coeffs + i, poly2->den, poly1->coeffs + i);
    for (i = min; i < poly2->length; i++)
    {
        fmpz_mul(res->coeffs + i, poly1->den, poly2->coeffs + i);
        fmpz_neg(res->coeffs + i, res->coeffs + i);
    }
    fmpz_mul(res->den, poly1->den, poly2->den);
    fmpz_clear(x);
    fmpz_clear(y);
    
    _fmpq_poly_set_length(res, max);
    fmpq_poly_canonicalise(res);
}

void fmpq_poly_sub(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
{
    ulong max = FLINT_MAX(poly1->length, poly2->length);
    fmpq_poly_fit_length(res, max);
    
    if (poly1 == poly2)
    {
        fmpq_poly_set_si(res, 0);
        return;
    }
    
    if (res == poly1)
        _fmpq_poly_sub_in_place(res->coeffs, res->den, res->length, 
                                poly2->coeffs, poly2->den, poly2->length);
    else if (res == poly2)
    {
        _fmpq_poly_sub_in_place(res->coeffs, res->den, res->length, 
                                poly1->coeffs, poly1->den, poly1->length);
        _fmpz_vec_neg(res->coeffs, res->coeffs, res->length);
    }
    else
        _fmpq_poly_sub(res->coeffs, res->den, 
                       poly1->coeffs, poly1->den, poly1->length, 
                       poly2->coeffs, poly2->den, poly2->length);
    
    _fmpq_poly_set_length(res, max);
    _fmpq_poly_normalise(res);
}

