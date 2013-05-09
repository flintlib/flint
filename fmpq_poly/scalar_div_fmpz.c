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

void _fmpq_poly_scalar_div_fmpz(fmpz * rpoly, fmpz_t rden, const fmpz * poly, 
                                const fmpz_t den, len_t len, const fmpz_t c)
{
    if (*c == 1L)
    {
        if (rpoly != poly)
        {
            _fmpz_vec_set(rpoly, poly, len);
            fmpz_set(rden, den);
        }
    }
    else if (*c == -1L)
    {
        _fmpz_vec_neg(rpoly, poly, len);
        fmpz_set(rden, den);
    }
    else
    {
        fmpz_t d;
        fmpz_init(d);
        _fmpz_vec_content(d, poly, len);
        fmpz_gcd(d, d, c);
        
        if (fmpz_sgn(c) < 0)
            fmpz_neg(d, d);
        _fmpz_vec_scalar_divexact_fmpz(rpoly, poly, len, d);
        fmpz_divexact(d, c, d);
        fmpz_mul(rden, den, d);
        
        fmpz_clear(d);
    }
}

void fmpq_poly_scalar_div_fmpz(fmpq_poly_t rop, const fmpq_poly_t op, const fmpz_t c)
{
    if (*c == 0L)
    {
        printf("Exception (fmpq_poly_scalar_div_fmpz). Division by zero.\n");
        abort();
    }
    
    if (fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }
    
    fmpq_poly_fit_length(rop, op->length);
    _fmpq_poly_set_length(rop, op->length);
    
    _fmpq_poly_scalar_div_fmpz(rop->coeffs, rop->den, 
                               op->coeffs, op->den, op->length, c);
}

