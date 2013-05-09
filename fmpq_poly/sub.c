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
#include "fmpq_poly.h"

void _fmpq_poly_sub(fmpz * rpoly, fmpz_t rden, 
                    const fmpz * poly1, const fmpz_t den1, len_t len1, 
                    const fmpz * poly2, const fmpz_t den2, len_t len2)
{
    len_t max = FLINT_MAX(len1, len2);
    len_t min = FLINT_MIN(len1, len2);
    
    fmpz_t d;
    fmpz_init(d);
    fmpz_one(d);
    if (*den1 != 1L && *den2 != 1L)
        fmpz_gcd(d, den1, den2);
    
    if (*d == 1L)
    {
        _fmpz_vec_scalar_mul_fmpz(rpoly, poly1, len1, den2);
        _fmpz_vec_scalar_submul_fmpz(rpoly, poly2, min, den1);
        if (len1 < len2)
        {
            _fmpz_vec_scalar_mul_fmpz(rpoly + min, poly2 + min, max - min, den1);
            _fmpz_vec_neg(rpoly + min, rpoly + min, max - min);
        }
        fmpz_mul(rden, den1, den2);
    }
    else
    {
        fmpz_t den11;
        fmpz_t den22;
        fmpz_init(den11);
        fmpz_init(den22);
        fmpz_divexact(den11, den1, d);
        fmpz_divexact(den22, den2, d);
        
        _fmpz_vec_scalar_mul_fmpz(rpoly, poly1, len1, den22);
        _fmpz_vec_scalar_submul_fmpz(rpoly, poly2, len2, den11);
        if (len1 < len2)
        {
            _fmpz_vec_scalar_mul_fmpz(rpoly + min, poly2 + min, max - min, den11);
            _fmpz_vec_neg(rpoly + min, rpoly + min, max - min);
        }
        
        if (_fmpz_vec_is_zero(rpoly, max))
            fmpz_one(rden);
        else
        {
            fmpz_t e;
            fmpz_init(e);
            _fmpz_vec_content(e, rpoly, max);
            if (*e != 1L)
                fmpz_gcd(e, e, d);
            
            if (*e == 1L)
                fmpz_mul(rden, den1, den22);
            else
            {
                _fmpz_vec_scalar_divexact_fmpz(rpoly, rpoly, max, e);
                fmpz_divexact(den11, den1, e);
                fmpz_mul(rden, den11, den22);
            }
            fmpz_clear(e);
        }
        fmpz_clear(den11);
        fmpz_clear(den22);
    }
    fmpz_clear(d);
}

void fmpq_poly_sub(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
{
    len_t len1, len2, max;

    if (poly1 == poly2)
    {
        fmpq_poly_zero(res);
        return;
    }
    
    len1 = poly1->length;
    len2 = poly2->length;
    max  = FLINT_MAX(poly1->length, poly2->length);
    fmpq_poly_fit_length(res, max);
    
    if (res != poly2)
        _fmpq_poly_sub(res->coeffs, res->den, 
                       poly1->coeffs, poly1->den, len1, 
                       poly2->coeffs, poly2->den, len2);
    else
    {
        _fmpq_poly_sub(res->coeffs, res->den, 
                       poly2->coeffs, poly2->den, len2, 
                       poly1->coeffs, poly1->den, len1);
        _fmpz_vec_neg(res->coeffs, res->coeffs, max);
    }
    
    _fmpq_poly_set_length(res, max);
    _fmpq_poly_normalise(res);
}

