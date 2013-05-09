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
#include "fmpq_poly.h"

void _fmpq_poly_scalar_div_fmpq(fmpz * rpoly, fmpz_t rden, 
                                const fmpz * poly, const fmpz_t den, len_t len, 
                                const fmpz_t r, const fmpz_t s)
{
    fmpz_t gcd1;  /* GCD( poly, r ) */
    fmpz_t gcd2;  /* GCD( s, den )  */
    fmpz_init(gcd1);
    fmpz_init(gcd2);
    fmpz_one(gcd1);
    fmpz_one(gcd2);
    if (*r != 1L)
    {
        _fmpz_vec_content(gcd1, poly, len);
        if (*gcd1 != 1L)
            fmpz_gcd(gcd1, gcd1, r);
    }
    if (*den != 1L && *s != 1L)
        fmpz_gcd(gcd2, s, den);
    
    if (*gcd1 == 1L)
    {
        if (*gcd2 == 1L)
        {
            _fmpz_vec_scalar_mul_fmpz(rpoly, poly, len, s);
            fmpz_mul(rden, den, r);
        }
        else
        {
            fmpz_t s2;
            fmpz_init(s2);
            fmpz_divexact(s2, s, gcd2);
            _fmpz_vec_scalar_mul_fmpz(rpoly, poly, len, s2);
            fmpz_divexact(rden, den, gcd2);
            fmpz_mul(rden, rden, r);
            fmpz_clear(s2);
        }
    }
    else
    {
        fmpz_t r2;
        fmpz_init(r2);
        fmpz_divexact(r2, r, gcd1);
        if (*gcd2 == 1L)
        {
            _fmpz_vec_scalar_divexact_fmpz(rpoly, poly, len, gcd1);
            _fmpz_vec_scalar_mul_fmpz(rpoly, rpoly, len, s);
            fmpz_mul(rden, den, r2);
        }
        else
        {
            fmpz_t s2;
            fmpz_init(s2);
            fmpz_divexact(s2, s, gcd2);
            _fmpz_vec_scalar_divexact_fmpz(rpoly, poly, len, gcd1);
            _fmpz_vec_scalar_mul_fmpz(rpoly, rpoly, len, s2);
            fmpz_divexact(rden, den, gcd2);
            fmpz_mul(rden, rden, r2);
            fmpz_clear(s2);
        }
        fmpz_clear(r2);
    }
    
    if (_fmpz_vec_is_zero(rpoly, len))
        fmpz_one(rden);
    if (fmpz_sgn(rden) < 0)
    {
        _fmpz_vec_neg(rpoly, rpoly, len);
        fmpz_neg(rden, rden);
    }
    
    fmpz_clear(gcd1);
    fmpz_clear(gcd2);
}

void fmpq_poly_scalar_div_fmpq(fmpq_poly_t rop, const fmpq_poly_t op, const fmpq_t c)
{
    if (fmpq_is_zero(c))
    {
        printf("Exception (fmpq_poly_scalar_div_fmpq). Division by zero.\n");
        abort();
    }

    if (fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
    }
    else
    {
        fmpq_poly_fit_length(rop, op->length);
        _fmpq_poly_set_length(rop, op->length);
        
        _fmpq_poly_scalar_div_fmpq(rop->coeffs, rop->den, 
                                   op->coeffs, op->den, op->length, 
                                   fmpq_numref(c), fmpq_denref(c));
    }
}

