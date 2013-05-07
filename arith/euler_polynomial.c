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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "arith.h"

void arith_euler_polynomial(fmpq_poly_t poly, ulong n)
{
    fmpz_t t;
    long k;

    if (n == 0)
    {
        fmpq_poly_set_ui(poly, 1UL);
        return;
    }

    arith_bernoulli_polynomial(poly, n + 1);

    fmpz_init(t);
    fmpz_set_si(t, -2L);
    for (k = n; k >= 0; k--)
    {
        fmpz_mul(poly->coeffs + k, poly->coeffs + k, t);
        fmpz_mul_ui(t, t, 2UL);
        fmpz_sub_ui(t, t, 2UL);
    }
    fmpz_zero(poly->coeffs + n + 1);
    fmpz_mul_ui(poly->den, poly->den, n + 1);
    fmpq_poly_canonicalise(poly);
    fmpz_clear(t);
}
