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

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"


void
_fmpz_poly_interpolate_fmpz_vec(fmpz * poly, fmpz_t den,
                                    const fmpz * xs, const fmpz * ys, long n)
{
    fmpz *P, *Q, *w, *linear;
    fmpz_t t;

    long i, j;

    /* Constant */
    if (n == 1)
    {
        fmpz_set(poly, ys);
        fmpz_set_ui(den, 1UL);
        return;
    }

    /* Linear */
    if (n == 2)
    {
        fmpz_sub(den, xs, xs + 1);
        fmpz_sub(poly + 1, ys, ys + 1);
        fmpz_mul(poly, xs, ys + 1);
        fmpz_submul(poly, xs + 1, ys);
        return;
    }

    fmpz_init(t);

    P = _fmpz_vec_init(n + 1);
    Q = _fmpz_vec_init(n);
    w = _fmpz_vec_init(n);
    linear = _fmpz_vec_init(2);

    /* P = (x-x[0])*(x-x[1])*...*(x-x[n-1]) */
    _fmpz_poly_product_roots_fmpz_vec(P, xs, n);

    /* Weights */
    for (i = 0; i < n; i++)
    {
        fmpz_set_ui(w + i, 1UL);
        for (j = 0; j < n; j++)
        {
            if (i != j)
            {
                fmpz_sub(t, xs + i, xs + j);
                fmpz_mul(w + i, w + i, t);
            }
        }
    }

    _fmpz_vec_zero(poly, n);
    _fmpz_vec_lcm(den, w, n);

    for (i = 0; i < n; i++)
    {
        /* Q = P / (x - x[i]) */
        fmpz_set_ui(linear + 1, 1UL);
        fmpz_neg(linear, xs + i);
        _fmpz_poly_div(Q, P, n + 1, linear, 2);

        /* result += Q * weight(i) */
        fmpz_divexact(t, den, w + i);
        fmpz_mul(t, t, ys + i);
        _fmpz_vec_scalar_addmul_fmpz(poly, Q, n, t);
    }

    _fmpz_vec_clear(P, n + 1);
    _fmpz_vec_clear(Q, n);
    _fmpz_vec_clear(w, n);
    _fmpz_vec_clear(linear, 2);
    fmpz_clear(t);
}

void
fmpz_poly_interpolate_fmpz_vec(fmpz_poly_t poly,
                                    const fmpz * xs, const fmpz * ys, long n)
{
    if (n == 0)
    {
        fmpz_poly_zero(poly);
        return;
    }
    else if (n == 1)
    {
        fmpz_poly_set_fmpz(poly, ys);
        return;
    }
    else
    {
        fmpz_t den;
        fmpz_init(den);
        fmpz_poly_fit_length(poly, n);
        _fmpz_poly_interpolate_fmpz_vec(poly->coeffs, den, xs, ys, n);
        _fmpz_vec_scalar_divexact_fmpz(poly->coeffs, poly->coeffs, n, den);
        _fmpz_poly_set_length(poly, n);
        _fmpz_poly_normalise(poly);
        fmpz_clear(den);
    }
}
