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

#include <stdlib.h>
#include <mpir.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void fmpq_poly_set_coeff_mpq(fmpq_poly_t poly, long n, const mpq_t x)
{
    long len = poly->length;
    const int replace = (n < len && !fmpz_is_zero(poly->coeffs + n));

    if (!replace && mpq_sgn(x) == 0)
        return;

    if (n + 1 > len)
    {
        fmpq_poly_fit_length(poly, n + 1);
        _fmpq_poly_set_length(poly, n + 1);
        mpn_zero((mp_ptr) poly->coeffs + len, (n + 1) - len);
        len = n + 1;
    }

    if (replace)
    {
        fmpz_t c;
        fmpz_t num, den;

        fmpz_init(c);
        fmpz_init(num);
        fmpz_init(den);

        fmpz_set_mpz(num, mpq_numref(x));
        fmpz_set_mpz(den, mpq_denref(x));

        fmpz_zero(poly->coeffs + n);
        _fmpz_poly_content(c, poly->coeffs, len);
        _fmpz_vec_scalar_mul_fmpz(poly->coeffs, poly->coeffs, len, den);
        fmpz_mul(c, c, den);

        fmpz_mul(poly->coeffs + n, num, poly->den);
        fmpz_gcd(c, c, poly->coeffs + n);
        fmpz_mul(poly->den, poly->den, den);

        if (*c != 1L)
            fmpz_gcd(c, c, poly->den);
        if (*c != 1L)
        {
            _fmpz_vec_scalar_divexact_fmpz(poly->coeffs, poly->coeffs, len, c);
            fmpz_divexact(poly->den, poly->den, c);
        }

        _fmpq_poly_normalise(poly);
        fmpz_clear(c);
        fmpz_clear(num);
        fmpz_clear(den);
    }
    else
    {
        fmpz_t d, t;

        fmpz_init(d);
        fmpz_init(t);
        fmpz_set_mpz(t, mpq_denref(x));

        fmpz_gcd(d, poly->den, t);
        fmpz_divexact(t, t, d);

        _fmpz_vec_scalar_mul_fmpz(poly->coeffs, poly->coeffs, len, t);
        fmpz_set_mpz(poly->coeffs + n, mpq_numref(x));
        fmpz_mul(poly->coeffs + n, poly->coeffs + n, poly->den);
        fmpz_divexact(poly->coeffs + n, poly->coeffs + n, d);

        fmpz_mul(poly->den, poly->den, t);

        fmpz_clear(d);
        fmpz_clear(t);
    }
}

