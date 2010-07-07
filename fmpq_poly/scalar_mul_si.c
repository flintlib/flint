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

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void fmpq_poly_scalar_mul_si(fmpq_poly_t rop, const fmpq_poly_t op, long c)
{
    fmpq_poly_fit_length(rop, op->length);
    
    if (*op->den == 1L)
    {
        _fmpz_vec_scalar_mul_si(rop->coeffs, op->coeffs, op->length, c);
        fmpz_set_si(rop->den, 1);
    }
    else
    {
        fmpz_t d, fc;
        fmpz_init(fc);
        fmpz_init(d);
        fmpz_set_si(fc, c);
        fmpz_gcd(d, op->den, fc);
        if (*d == 1L)
        {
            _fmpz_vec_scalar_mul_si(rop->coeffs, op->coeffs, op->length, c);
            fmpz_set(rop->den, op->den);
        }
        else
        {
            fmpz_divexact(fc, fc, d);
            _fmpz_vec_scalar_mul_fmpz(rop->coeffs, op->coeffs, op->length, fc);
            fmpz_divexact(rop->den, op->den, d);
        }
        
        fmpz_clear(fc);
        fmpz_clear(d);
    }
    
    _fmpq_poly_set_length(rop, op->length);
    _fmpq_poly_normalise(rop);
}

