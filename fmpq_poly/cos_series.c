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
                       const fmpz * h, const fmpz_t hden, slong hlen, slong n)
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
    fmpz_mul_ui(uden, hden, UWORD(2));
    _fmpq_poly_tan_series(t, tden, h, uden, hlen, n);
    _fmpq_poly_mullow(u, uden, t, tden, n, t, tden, n, n);

    _fmpq_poly_canonicalise(u, uden, n);
    _fmpz_vec_neg(t, u, n);
    fmpz_set(tden, uden);
    fmpz_set(t, tden);
    fmpz_set(u, uden);

    _fmpq_poly_div_series(g, gden, t, tden, n, u, uden, n, n);
    _fmpq_poly_canonicalise(g, gden, n);

    _fmpz_vec_clear(t, n);
    _fmpz_vec_clear(u, n);
    fmpz_clear(tden);
    fmpz_clear(uden);
}

void fmpq_poly_cos_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
{
    if (n == 0)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (poly->length == 0 || n == 1)
    {
        fmpq_poly_one(res);
        return;
    }

    if (!fmpz_is_zero(poly->coeffs))
    {
        flint_printf("Exception (fmpq_poly_cos_series). Constant term != 0.\n");
        abort();
    }

    if (res != poly)
    {
        fmpq_poly_fit_length(res, n);
        _fmpq_poly_cos_series(res->coeffs, res->den,
            poly->coeffs, poly->den, poly->length, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_cos_series(t->coeffs, t->den,
            poly->coeffs, poly->den, poly->length, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);
}
