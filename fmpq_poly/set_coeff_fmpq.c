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

    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2010 William Hart

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void fmpq_poly_set_coeff_fmpq(fmpq_poly_t poly, long n, const fmpq_t x)
{
    long len = poly->length;
    const int replace = (n < len && !fmpz_is_zero(poly->coeffs + n));

    if (!replace && fmpq_is_zero(x))
        return;

    if (n + 1 > len)
    {
        fmpq_poly_fit_length(poly, n + 1);
        _fmpq_poly_set_length(poly, n + 1);
        flint_mpn_zero((mp_ptr) poly->coeffs + len, (n + 1) - len);
        len = n + 1;
    }

    if (replace)
    {
        fmpz_t c;

        fmpz_init(c);

        fmpz_zero(poly->coeffs + n);
        _fmpz_poly_content(c, poly->coeffs, len);
        _fmpz_vec_scalar_mul_fmpz(poly->coeffs, poly->coeffs, len, fmpq_denref(x));
        fmpz_mul(c, c, fmpq_denref(x));

        fmpz_mul(poly->coeffs + n, fmpq_numref(x), poly->den);
        fmpz_gcd(c, c, poly->coeffs + n);
        fmpz_mul(poly->den, poly->den, fmpq_denref(x));

        if (!fmpz_is_one(c))
            fmpz_gcd(c, c, poly->den);
        if (!fmpz_is_one(c))
        {
            _fmpz_vec_scalar_divexact_fmpz(poly->coeffs, poly->coeffs, len, c);
            fmpz_divexact(poly->den, poly->den, c);
        }

        _fmpq_poly_normalise(poly);
        fmpz_clear(c);
    }
    else
    {
        fmpz_t d, t;

        fmpz_init(d);
        fmpz_init(t);

        fmpz_gcd(d, poly->den, fmpq_denref(x));
        fmpz_divexact(t, fmpq_denref(x), d);

        _fmpz_vec_scalar_mul_fmpz(poly->coeffs, poly->coeffs, len, t);
        fmpz_set(poly->coeffs + n, fmpq_numref(x));
        fmpz_mul(poly->coeffs + n, poly->coeffs + n, poly->den);
        fmpz_divexact(poly->coeffs + n, poly->coeffs + n, d);

        fmpz_mul(poly->den, poly->den, t);

        fmpz_clear(d);
        fmpz_clear(t);
    }
}

