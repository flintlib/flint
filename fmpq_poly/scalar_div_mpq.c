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
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void _fmpq_poly_scalar_div_mpq(fmpz * rpoly, fmpz_t rden, 
                               const fmpz * poly, const fmpz_t den, long len, 
                               const fmpz_t r, const fmpz_t s)
{
    fmpz_t gcd1;  /* GCD( poly, r ) */
    fmpz_t gcd2;  /* GCD( s, den )  */
    fmpz_init(gcd1);
    fmpz_init(gcd2);
    fmpz_set_ui(gcd1, 1);
    fmpz_set_ui(gcd2, 1);
    if (*r != 1L)
    {
        _fmpz_vec_content(gcd1, poly, len);
        if (*gcd1 != 1L)
            fmpz_gcd(gcd1, gcd1, r);
    }
    if (*den != 1L && *s != 1L)
        fmpz_gcd(gcd2, gcd2, den);
    
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
        fmpz_set_ui(rden, 1);
    if (fmpz_sgn(rden) < 0)
    {
        _fmpz_vec_neg(rpoly, rpoly, len);
        fmpz_neg(rden, rden);
    }
    
    fmpz_clear(gcd1);
    fmpz_clear(gcd2);
}

void fmpq_poly_scalar_div_mpq(fmpq_poly_t rop, const fmpq_poly_t op, const mpq_t c)
{
    fmpz_t r, s;

    if (mpq_sgn(c) == 0)
    {
        printf("Exception: division by zero in fmpq_poly_scalar_div_mpq\n");
        abort();
    }
    
    if (fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }
    
    fmpq_poly_fit_length(rop, op->length);
    _fmpq_poly_set_length(rop, op->length);
    
    fmpz_init(r);
    fmpz_init(s);
    fmpz_set_mpz(r, mpq_numref(c));
    fmpz_set_mpz(s, mpq_denref(c));
    
    _fmpq_poly_scalar_div_mpq(rop->coeffs, rop->den, 
                              op->coeffs, op->den, op->length, r, s);
    
    fmpz_clear(r);
    fmpz_clear(s);
}

