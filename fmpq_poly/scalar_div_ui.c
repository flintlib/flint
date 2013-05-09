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

void _fmpq_poly_scalar_div_ui(fmpz * rpoly, fmpz_t rden, const fmpz * poly, 
                              const fmpz_t den, len_t len, ulong c)
{
    if (c == 1UL)
    {
        if (rpoly != poly)
            _fmpz_vec_set(rpoly, poly, len);
        fmpz_set(rden, den);
    }
    else
    {
        fmpz_t d, fc;
        ulong ud;
        fmpz_init(d);
        fmpz_init(fc);
        _fmpz_vec_content(d, poly, len);
        fmpz_set_ui(fc, c);
        fmpz_gcd(d, d, fc);
        ud = fmpz_get_ui(d);  /* gcd of d and c fits into a ulong */
        
        _fmpz_vec_scalar_divexact_ui(rpoly, poly, len, ud);
        fmpz_mul_ui(rden, den, c / ud);
        
        fmpz_clear(d);
        fmpz_clear(fc);
    }
}

void fmpq_poly_scalar_div_ui(fmpq_poly_t rop, const fmpq_poly_t op, ulong c)
{
    if (c == 0UL)
    {
        printf("Exception (fmpq_poly_scalar_div_ui). Division by zero.\n");
        abort();
    }
    
    if (fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }
    
    fmpq_poly_fit_length(rop, op->length);
    _fmpq_poly_set_length(rop, op->length);
    
    _fmpq_poly_scalar_div_ui(rop->coeffs, rop->den, 
                             op->coeffs, op->den, op->length, c);
}

