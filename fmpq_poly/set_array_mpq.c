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

    Copyright (C) 2010, 2011 Sebastian Pancratz

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void
_fmpq_poly_set_array_mpq(fmpz * poly, fmpz_t den, const mpq_t * a, len_t n)
{
    len_t i;
    mpz_t d, t;

    mpz_init_set_ui(d, 1);
    mpz_init(t);
    for (i = 0; i < n; i++)
    {
        mpz_lcm(d, d, mpq_denref(a[i]));
    }

    for (i = 0; i < n; i++)
    {
        __mpz_struct *ptr = _fmpz_promote(poly + i);

        mpz_divexact(t, d, mpq_denref(a[i]));
        mpz_mul(ptr, mpq_numref(a[i]), t);
        _fmpz_demote_val(poly + i);
    }

    fmpz_set_mpz(den, d);
    mpz_clear(d);
    mpz_clear(t);
}

void fmpq_poly_set_array_mpq(fmpq_poly_t poly, const mpq_t * a, len_t n)
{
    if (n == 0)
    {
        fmpq_poly_zero(poly);
    }
    else
    {
        fmpq_poly_fit_length(poly, n);
        _fmpq_poly_set_array_mpq(poly->coeffs, poly->den, a, n);
        _fmpq_poly_set_length(poly, n);
        _fmpq_poly_normalise(poly);
    }
}
