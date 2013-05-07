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
#include "fmpz_poly.h"
#include "fmpq_poly.h"


static void
_set_vec(fmpz * rnum, fmpz_t den,
                const fmpz * xnum, const fmpz * xden, long len)
{
    long j;
    fmpz_t t;
    fmpz_init(t);
    fmpz_one(den);

    for (j = 0; j < len; j++)
        fmpz_lcm(den, den, xden + j);

    for (j = 0; j < len; j++)
    {
        fmpz_divexact(t, den, xden + j);
        fmpz_mul(rnum + j, xnum + j, t);
    }

    fmpz_clear(t);
}

void
_fmpq_poly_revert_series_lagrange(fmpz * Qinv, fmpz_t den,
                            const fmpz * Q, const fmpz_t Qden, long n)
{
    long i;
    fmpz *R, *S, *T, *dens, *tmp;
    fmpz_t Rden, Sden, Tden;

    if (fmpz_is_one(Qden) && (n > 1) && fmpz_is_pm1(Q + 1))
    {
        _fmpz_poly_revert_series(Qinv, Q, n);
        fmpz_one(den);
    }
    else if (n <= 2)
    {
        fmpz_zero(Qinv);
        if (n == 2)
        {
            fmpz_set(Qinv + 1, Qden);
            fmpz_set(den, Q + 1);
            _fmpq_poly_canonicalise(Qinv, den, 2);
        }
    }
    else
    {
        dens = _fmpz_vec_init(n);
        R = _fmpz_vec_init(n - 1);
        S = _fmpz_vec_init(n - 1);
        T = _fmpz_vec_init(n - 1);
        fmpz_init(Rden);
        fmpz_init(Sden);
        fmpz_init(Tden);

        fmpz_zero(Qinv);
        fmpz_one(dens);
        fmpz_set(Qinv + 1, Qden);
        fmpz_set(dens + 1, Q + 1);

        _fmpq_poly_inv_series(R, Rden, Q + 1, Qden, n - 1);
        _fmpq_poly_canonicalise(R, Rden, n - 1);

        _fmpz_vec_set(S, R, n - 1);
        fmpz_set(Sden, Rden);

        for (i = 2; i < n; i++)
        {
            _fmpq_poly_mullow(T, Tden, S, Sden, n - 1, R, Rden, n - 1, n - 1);
            _fmpq_poly_canonicalise(T, Tden, n - 1);
            fmpz_set(Qinv + i, T + i - 1);
            fmpz_mul_ui(dens + i, Tden, i);
            tmp = S; S = T; T = tmp;
            fmpz_swap(Sden, Tden);
        }

        _set_vec(Qinv, den, Qinv, dens, n);
        _fmpq_poly_canonicalise(Qinv, den, n);

        _fmpz_vec_clear(R, n - 1);
        _fmpz_vec_clear(S, n - 1);
        _fmpz_vec_clear(T, n - 1);
        _fmpz_vec_clear(dens, n);
        fmpz_clear(Rden);
        fmpz_clear(Sden);
        fmpz_clear(Tden);
    }
}

void
fmpq_poly_revert_series_lagrange(fmpq_poly_t res,
            const fmpq_poly_t poly, long n)
{
    fmpz *copy;
    int alloc;

    if (poly->length < 2 || !fmpz_is_zero(poly->coeffs)
                         || fmpz_is_zero(poly->coeffs + 1))
    {
        printf("Exception (fmpq_poly_revert_series_lagrange). Input must have \n"
               "zero constant term and nonzero coefficient of x^1.\n");
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
        _fmpq_poly_revert_series_lagrange(res->coeffs,
                res->den, copy, poly->den, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_revert_series_lagrange(t->coeffs,
                t->den, copy, poly->den, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);

    if (alloc)
        flint_free(copy);
}
