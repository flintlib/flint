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

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void _fmpq_poly_log_series(fmpz * g, fmpz_t gden, 
                       const fmpz * f, const fmpz_t fden, slong flen, slong n)
{
    fmpz * f_diff;
    fmpz * f_inv;
    fmpz_t f_diff_den;
    fmpz_t f_inv_den;

    flen = FLINT_MIN(flen, n);

    f_diff = _fmpz_vec_init(flen - 1);
    f_inv = _fmpz_vec_init(n);
    fmpz_init(f_diff_den);
    fmpz_init(f_inv_den);

    _fmpq_poly_derivative(f_diff, f_diff_den, f, fden, flen);
    _fmpq_poly_inv_series(f_inv, f_inv_den, f, fden, flen, n);
    _fmpq_poly_mullow(g, gden, f_inv, f_inv_den, n - 1,
        f_diff, f_diff_den, flen - 1, n - 1);
    _fmpq_poly_integral(g, gden, g, gden, n);

    _fmpz_vec_clear(f_diff, flen - 1);
    _fmpz_vec_clear(f_inv, n);
    fmpz_clear(f_diff_den);
    fmpz_clear(f_inv_den);
}

void
fmpq_poly_log_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)
{
    slong flen = f->length;

    if (flen < 1 || !fmpz_equal(f->coeffs, f->den))
    {
        flint_printf("Exception (fmpq_poly_log_series). Constant term != 1.\n");
        flint_abort();
    }

    if (flen == 1 || n < 2)
    {
        fmpq_poly_zero(res);
        return;
    }

    fmpq_poly_fit_length(res, n);
    _fmpq_poly_log_series(res->coeffs, res->den, f->coeffs, f->den, f->length, n);
    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);
}
