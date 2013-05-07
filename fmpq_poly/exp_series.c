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
_fmpq_poly_exp_series(fmpz * g, fmpz_t gden,
                        const fmpz * h, const fmpz_t hden, long n)
{
    long m;
    fmpz * t, * u;
    fmpz_t tden, uden;

    if (n < 2)
    {
        fmpz_one(g);
        fmpz_one(gden);
        return;
    }

    m = (n + 1) / 2;
    _fmpq_poly_exp_series(g, gden, h, hden, m);
    _fmpz_vec_zero(g + m, n - m);

    t = _fmpz_vec_init(n);
    u = _fmpz_vec_init(n);
    fmpz_init(tden);
    fmpz_init(uden);

    _fmpq_poly_log_series(t, tden, g, gden, n);
    _fmpq_poly_sub(t, tden, t, tden, n, h, hden, n);
    /* TODO: half of product is redundant! */
    _fmpq_poly_mullow(u, uden, g, gden, n, t, tden, n, n);
    _fmpq_poly_sub(g, gden, g, gden, n, u, uden, n);
    _fmpq_poly_canonicalise(g, gden, n);

    fmpz_clear(tden);
    fmpz_clear(uden);
    _fmpz_vec_clear(t, n);
    _fmpz_vec_clear(u, n);
}

void fmpq_poly_exp_series(fmpq_poly_t res, const fmpq_poly_t poly, long n)
{
    fmpz *copy;
    int alloc;

    if (poly->length == 0)
    {
        fmpq_poly_set_ui(res, 1UL);
        return;
    }

    if (!fmpz_is_zero(poly->coeffs))
    {
        printf("Exception (fmpq_poly_exp_series). Constant term != 0.\n");
        abort();
    }

    if (n < 2)
    {
        if (n == 0) fmpq_poly_zero(res);
        if (n == 1) fmpq_poly_set_ui(res, 1UL);
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
        _fmpq_poly_exp_series(res->coeffs, res->den, copy, poly->den, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_exp_series(t->coeffs, t->den, copy, poly->den, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);

    if (alloc)
        flint_free(copy);
}

