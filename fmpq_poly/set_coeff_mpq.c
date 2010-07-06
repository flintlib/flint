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
    int canonicalise;
    fmpz_t s, t;
    
    /* We only need to canonicalise at the end if we replace a zero          */
    /* coefficient.                                                          */
    canonicalise = (n < poly->length) && !(fmpz_is_zero(poly->coeffs + n));
    
    /* Ensure there is enough space, and insert zeroes between the           */
    /* end of poly and the new coefficient if needed                         */
    fmpq_poly_fit_length(poly, n + 1);
    if (n + 1 > poly->length)
    {
        mpn_zero(poly->coeffs + poly->length, n - poly->length);
        poly->length = n + 1;
    }
    
    fmpz_init(s);
    fmpz_init(t);
    fmpz_set_mpz(s, mpq_denref(x));
    fmpz_gcd(t, s, poly->den);
    
    if (fmpz_equal(s, t))
    {
        fmpz_set_mpz(s, mpq_numref(x));
        fmpz_divexact(t, poly->den, t);
        fmpz_mul(s, s, t);
        fmpz_set(poly->coeffs + n, s);
    }
    else
    {
        fmpz_divexact(s, s, t);
        _fmpz_vec_scalar_mul_fmpz(poly->coeffs, poly->coeffs, poly->length, s);
        
        fmpz_divexact(poly->coeffs + n, poly->den, t);
        fmpz_set_mpz(t, mpq_numref(x));
        fmpz_mul(poly->coeffs + n, poly->coeffs + n, t);
        
        fmpz_mul(poly->den, poly->den, s);
    }
    
    if (canonicalise)
        fmpq_poly_canonicalise(poly);
    
    fmpz_clear(s);
    fmpz_clear(t);
}

