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

#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void _fmpq_poly_scalar_mul_si(fmpz * rpoly, fmpz_t rden, 
                              const fmpz * poly, const fmpz_t den, long len, 
                              long c)
{
    fmpz_t gcd;  /* GCD( den, c ) */

    if (c == 0)
    {
        _fmpz_vec_zero(rpoly, len);
        fmpz_one(rden);
        return;
    }

    fmpz_init(gcd);
    fmpz_set_si(gcd, c);
    fmpz_gcd(gcd, gcd, den);
    if (*gcd == 1L)
    {
        _fmpz_vec_scalar_mul_si(rpoly, poly, len, c);
        fmpz_set(rden, den);
    }
    else
    {
        if (c > LONG_MIN || fmpz_cmp_ui(gcd, - (ulong) LONG_MIN))
        {
            long g = fmpz_get_si(gcd);

            _fmpz_vec_scalar_mul_si(rpoly, poly, len, c / g);
            fmpz_divexact_si(rden, den, g);
        }
        else
        {
            _fmpz_vec_neg(rpoly, poly, len);
            fmpz_divexact_ui(rden, den, - (ulong) LONG_MIN);
        }
    }
    fmpz_clear(gcd);
}

void fmpq_poly_scalar_mul_si(fmpq_poly_t rop, const fmpq_poly_t op, long c)
{
    if (c == 0 || fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }
    
    fmpq_poly_fit_length(rop, op->length);
    _fmpq_poly_set_length(rop, op->length);
    
    _fmpq_poly_scalar_mul_si(rop->coeffs, rop->den, 
                             op->coeffs, op->den, op->length, c);
}

