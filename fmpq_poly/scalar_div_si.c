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

void _fmpq_poly_scalar_div_si(fmpz * rpoly, fmpz_t rden, const fmpz * poly,
                              const fmpz_t den, long len, long c)
{
    if (c == 1)
    {
        if (rpoly != poly)
        {
            _fmpz_vec_set(rpoly, poly, len);
            fmpz_set(rden, den);
        }
    }
    else if (c == -1)
    {
        _fmpz_vec_neg(rpoly, poly, len);
        fmpz_set(rden, den);
    }
    else
    {
        fmpz_t d, f;

        fmpz_init(d);
        fmpz_init(f);
        
        fmpz_set_si(f, c);
        _fmpz_vec_content(d, poly, len);
        fmpz_gcd(d, d, f);

        if (c > 0)
        {
            _fmpz_vec_scalar_divexact_fmpz(rpoly, poly, len, d);
            fmpz_mul_si(rden, den, c / fmpz_get_si(d));
        }
        else
        {
            ulong q = (- (ulong) c) / fmpz_get_ui(d);

            fmpz_neg(d, d);
            _fmpz_vec_scalar_divexact_fmpz(rpoly, poly, len, d);
            fmpz_mul_ui(rden, den, q);
        }

        fmpz_clear(d);
        fmpz_clear(f);
    }
}

void fmpq_poly_scalar_div_si(fmpq_poly_t rop, const fmpq_poly_t op, long c)
{
    if (c == 0L)
    {
        printf("Exception (fmpq_poly_scalar_div_si). Division by zero.\n");
        abort();
    }
    
    if (fmpq_poly_is_zero(op))
    {
        fmpq_poly_zero(rop);
        return;
    }
    
    fmpq_poly_fit_length(rop, op->length);
    _fmpq_poly_set_length(rop, op->length);
    
    _fmpq_poly_scalar_div_si(rop->coeffs, rop->den, 
                             op->coeffs, op->den, op->length, c);
}

