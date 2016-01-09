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

  Copyright (C) 2016 Shivin Srivastava

 ******************************************************************************/

#include "arith.h"

void _arith_fibonacci_polynomial(fmpz * coeffs, ulong n)
{
	fmpz * r;
	int even;
	slong k;
	ulong L;

	L = n / 2;
	even = 1 - (n % 2);
	/* set the first two coefficients of poly depending parity of n */
	if (even)
	{
		fmpz_zero(coeffs);
		fmpz_one(coeffs + 1);
		fmpz_mul_ui(coeffs + 1, coeffs + 1, L);
	}

	else 
	{
		fmpz_one(coeffs);
		fmpz_zero(coeffs + 1);
	}

	r = coeffs + even;
	r += 2;

	/* calculate the coefficients of the polynomial*/
	for (k = 2 + even; k < n; k += 2)
	{
		fmpz_mul2_uiui(r, r - 2, L + k / 2, L + k / 2 - k + 1);
		fmpz_divexact2_uiui(r, r, k, k - 1);
		r += 2;
	}

	/* set the alternate coefficients to 0 again depending on the parity*/
	for (k = 1 + even; k < n; k += 2)
	{
		fmpz_zero(coeffs + k);
	}
}


void arith_fibonacci_polynomial(fmpz_poly_t poly, ulong n)
{
	if (n == 0)
	{
		fmpz_poly_set_ui(poly, UWORD(0));
		return;
	}

	if (n == 1)
	{
		fmpz_poly_set_ui(poly, UWORD(1));
	}

	else
	{
		fmpz_poly_fit_length(poly, n);
		_arith_fibonacci_polynomial(poly->coeffs, n);
	}
	_fmpz_poly_set_length(poly, n);
}

