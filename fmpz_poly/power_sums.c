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

    Copyright (C) 2016 Vincent Delecroix

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void
fmpz_poly_power_sums(fmpz_poly_t res, const fmpz_poly_t poly, slong n)
{
    slong len = poly->length;

    if (len == 0)
    {
        flint_printf
            ("Exception (fmpz_poly_power_sums). Zero polynomial.\n");
        abort();
    }
    else if (n <= 0 || len == 1)
    {
        fmpz_poly_zero(res);
    }
    else
    {
        size_t i = 0;
        while (fmpz_is_zero(poly->coeffs + i))
            i++;
        if (poly == res)
        {
            fmpz_poly_t t;
            fmpz_poly_init2(t, n);
            _fmpz_poly_power_sums_naive(t->coeffs, poly->coeffs + i,
                                           len - i, n);
            fmpz_poly_swap(res, t);
            fmpz_poly_clear(t);
        }
        else
        {
            fmpz_poly_fit_length(res, n);
            _fmpz_poly_power_sums_naive(res->coeffs, poly->coeffs + i,
                                           len - i, n);
        }
        _fmpz_poly_set_length(res, n);
        if (i)
            fmpz_set_si(res->coeffs, len - 1);
        _fmpz_poly_normalise(res);
    }
}
