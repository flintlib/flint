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

void _fmpq_poly_scalar_mul_ui(fmpz * rpoly, fmpz_t rden, 
                              const fmpz * poly, const fmpz_t den, len_t len, 
                              ulong c)
{
    fmpz_t gcd;  /* GCD( den, c ) */

    if (c == 0)
    {
        _fmpz_vec_zero(rpoly, len);
        fmpz_one(rden);
        return;
    }

    fmpz_init(gcd);
    fmpz_set_ui(gcd, c);
    fmpz_gcd(gcd, gcd, den);
    if (*gcd == 1L)
    {
        _fmpz_vec_scalar_mul_ui(rpoly, poly, len, c);
        fmpz_set(rden, den);
    }
    else
    {
        ulong gcd2 = fmpz_get_ui(gcd);
        ulong c2 = c / gcd2;
        _fmpz_vec_scalar_mul_ui(rpoly, poly, len, c2);
        fmpz_fdiv_q_ui(rden, den, gcd2);
    }
    fmpz_clear(gcd);
}

void fmpq_poly_scalar_mul_ui(fmpq_poly_t rop, const fmpq_poly_t op, ulong c)
{
    if (c == 0 || fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }
    
    fmpq_poly_fit_length(rop, op->length);
    _fmpq_poly_set_length(rop, op->length);
    
    _fmpq_poly_scalar_mul_ui(rop->coeffs, rop->den, 
                             op->coeffs, op->den, op->length, c);
}

