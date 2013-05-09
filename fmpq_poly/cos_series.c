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

void
_fmpq_poly_cos_series(fmpz * g, fmpz_t gden, 
                           const fmpz * h, const fmpz_t hden, len_t n)
{
    fmpz * t;
    fmpz * u;
    fmpz_t tden;
    fmpz_t uden;

    t = _fmpz_vec_init(n);
    u = _fmpz_vec_init(n);
    fmpz_init(tden);
    fmpz_init(uden);

    /* cos(x) = (1-tan(x/2)^2)/(1+tan(x/2)^2) */
    fmpz_mul_ui(uden, hden, 2UL);
    _fmpq_poly_tan_series(t, tden, h, uden, n);
    _fmpq_poly_mullow(u, uden, t, tden, n, t, tden, n, n);

    _fmpq_poly_canonicalise(u, uden, n);
    _fmpz_vec_neg(t, u, n);
    fmpz_set(tden, uden);
    fmpz_set(t, tden);
    fmpz_set(u, uden);

    _fmpq_poly_div_series(g, gden, t, tden, u, uden, n);
    _fmpq_poly_canonicalise(g, gden, n);

    _fmpz_vec_clear(t, n);
    _fmpz_vec_clear(u, n);
    fmpz_clear(tden);
    fmpz_clear(uden);
}

void
fmpq_poly_cos_series(fmpq_poly_t res, const fmpq_poly_t f, len_t n)
{
    fmpz * f_coeffs;
    len_t flen = f->length;

    if (f->length == 0)
    {
        fmpq_poly_set_ui(res, 1UL);
        return;
    }

    if (!fmpz_is_zero(f->coeffs))
    {
        printf("Exception (fmpq_poly_cos_series). Constant term != 0.");
        abort();
    }

    if (n < 2)
    {
        if (n == 0) fmpq_poly_zero(res);
        if (n == 1) fmpq_poly_set_ui(res, 1UL);
        return;
    }

    fmpq_poly_fit_length(res, n);

    if (flen < n)
    {
        f_coeffs = _fmpz_vec_init(n);
        _fmpz_vec_set(f_coeffs, f->coeffs, flen);
    }
    else
    {
        f_coeffs = f->coeffs;
    }

    _fmpq_poly_cos_series(res->coeffs, res->den, f_coeffs, f->den, n);

    if (flen < n)
    {
        _fmpz_vec_clear(f_coeffs, n);
    }

    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);
}
