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

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void
_fmpq_poly_tan_series(fmpz * g, fmpz_t gden,
                        const fmpz * h, const fmpz_t hden, long n)
{
    long m;
    fmpz * t, * u, * v;
    fmpz_t tden, uden, vden;

    if (n <= 3)
    {
        if (n >= 1) fmpz_set(g + 0, h + 0);
        if (n >= 2) fmpz_set(g + 1, h + 1);
        if (n == 3) fmpz_set(g + 2, h + 2);
        fmpz_set(gden, hden);
        _fmpq_poly_canonicalise(g, gden, n);
        return;
    }

    m = (n + 1) / 2;

    _fmpq_poly_tan_series(g, gden, h, hden, m);
    _fmpz_vec_zero(g + m, n - m);

    t = _fmpz_vec_init(n);
    u = _fmpz_vec_init(n);
    v = _fmpz_vec_init(n);
    fmpz_init(tden);
    fmpz_init(uden);
    fmpz_init(vden);

    _fmpq_poly_mul(u, uden, g, gden, m, g, gden, m);
    fmpz_set(u, uden); /* u += 1 */
    if (2*m - 1 < n)
        fmpz_zero(u + n - 1);

    _fmpq_poly_atan_series(t, tden, g, gden, n);
    _fmpq_poly_sub(t, tden, t, tden, n, h, hden, n);
    _fmpq_poly_mullow(v + m, vden, u, uden, n, t + m, tden, n - m, n - m);
    _fmpq_poly_sub(g, gden, g, gden, m, v, vden, n);
    _fmpq_poly_canonicalise(g, gden, n);

    fmpz_clear(tden);
    fmpz_clear(uden);
    fmpz_clear(vden);
    _fmpz_vec_clear(t, n);
    _fmpz_vec_clear(u, n);
    _fmpz_vec_clear(v, n);
}

void fmpq_poly_tan_series(fmpq_poly_t res, const fmpq_poly_t poly, long n)
{
    fmpz *copy;
    int alloc;

    if (poly->length == 0)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (!fmpz_is_zero(poly->coeffs))
    {
        printf("Exception (fmpq_poly_tan_series). Constant term != 0.\n");
        abort();
    }

    if (n < 2)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (poly->length >= n)
    {
        copy = poly->coeffs;
        alloc = 0;
    }
    else
    {
        long i;
        copy = (fmpz *) flint_malloc(n * sizeof(fmpz));
        for (i = 0; i < poly->length; i++)
            copy[i] = poly->coeffs[i];
        for ( ; i < n; i++)
            copy[i] = 0;
        alloc = 1;
    }

    if (res != poly)
    {
        fmpq_poly_fit_length(res, n);
        _fmpq_poly_tan_series(res->coeffs, res->den, copy, poly->den, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_tan_series(t->coeffs, t->den, copy, poly->den, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);

    if (alloc)
        flint_free(copy);
}

