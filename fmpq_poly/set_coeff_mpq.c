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
    Copyright (C) 2010 William Hart

******************************************************************************/

#include <mpir.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void fmpq_poly_set_coeff_mpq(fmpq_poly_t poly, ulong n, const mpq_t x)
{
    ulong len = poly->length;
    int replace = (n < len && *(poly->coeffs + n) != 0L);
    
    if (!replace && mpq_sgn(x) == 0)
        return;
    
    if (n + 1UL > len)
    {
        fmpq_poly_fit_length(poly, n + 1UL); 
        mpn_zero(poly->coeffs + len, n - len);
        len = n + 1UL;
        poly->length = len;
    }
    
    if (replace)
    {
        fmpz_t c;
        fmpz_t t;
        fmpz_init(c);
        fmpz_init(t);
        
        fmpz_set_mpz(poly->coeffs + n, mpq_numref(x));
        _fmpz_vec_content(c, poly->coeffs, len);
        if (*c != 1L)
            fmpz_gcd(c, c, poly->den);
        if (*c != 1L)
        {
            _fmpz_vec_scalar_divexact(poly->coeffs, poly->coeffs, len, c);
            fmpz_divexact(poly->den, poly->den, c);
        }
        fmpz_set_mpz(t, mpq_denref(x));
        _fmpz_vec_scalar_mul_fmpz(poly->coeffs, poly->coeffs, len, t);
        fmpz_set_mpz(poly->coeffs + n, mpq_numref(x));
        fmpz_mul(poly->coeffs + n, poly->coeffs + n, poly->den);
        fmpz_mul(poly->den, t, poly->den);
        
        _fmpq_poly_normalise(poly);
        fmpz_clear(c);
        fmpz_clear(t);
    }
    else
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_set_mpz(t, mpq_denref(x));
        _fmpz_vec_scalar_mul_fmpz(poly->coeffs, poly->coeffs, len, t);
        fmpz_set_mpz(poly->coeffs + n, mpq_numref(x));
        fmpz_mul(poly->coeffs + n, poly->coeffs + n, poly->den);
        fmpz_mul(poly->den, t, poly->den);
        fmpz_clear(t);
    }
}

