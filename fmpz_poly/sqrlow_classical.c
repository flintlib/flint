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

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

/*
    Assumes len > 0 and 0 < n <= 2 * len - 1.
 */
void _fmpz_poly_sqrlow_classical(fmpz *rop, const fmpz *op, len_t len, len_t n)
{
    if (len == 1 || n == 1)  /* Special case */
    {
        fmpz_mul(rop, op, op);
    }
    else   /* Ordinary case */
    {
        len_t i;

        _fmpz_vec_scalar_mul_fmpz(rop, op, FLINT_MIN(len, n), op);

        _fmpz_vec_scalar_mul_fmpz(rop + len, op + 1, n - len, op + len - 1);

        for (i = 1; i < len - 1; i++)
            _fmpz_vec_scalar_addmul_fmpz(rop + i + 1, 
                op + 1, FLINT_MIN(i - 1, n - (i + 1)), op + i);

        for (i = 1; i < FLINT_MIN(2 * len - 2, n); i++)
            fmpz_mul_ui(rop + i, rop + i, 2);

        for (i = 1; i < FLINT_MIN(len - 1, (n + 1) / 2); i++)
            fmpz_addmul(rop + 2 * i, op + i, op + i);
    }
}

void
fmpz_poly_sqrlow_classical(fmpz_poly_t res, const fmpz_poly_t poly, len_t n)
{
    len_t len = poly->length;

    if (len == 0 || n == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    n = FLINT_MIN(2 * len - 1, n);

    if (res == poly)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        _fmpz_poly_sqrlow_classical(t->coeffs, poly->coeffs, len, n);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
    else
    {
        fmpz_poly_fit_length(res, n);
        _fmpz_poly_sqrlow_classical(res->coeffs, poly->coeffs, len, n);
    }

    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);
}

