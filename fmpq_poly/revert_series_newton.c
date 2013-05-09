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


#define FLINT_REVERSE_NEWTON_CUTOFF 4

void
_fmpq_poly_revert_series_newton(fmpz * Qinv, fmpz_t den,
        const fmpz * Q, const fmpz_t Qden, len_t n)
{
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
        len_t *a, i, k;
        fmpz *T, *U, *V;
        fmpz_t Tden, Uden, Vden;

        T = _fmpz_vec_init(n);
        U = _fmpz_vec_init(n);
        V = _fmpz_vec_init(n);
        fmpz_init(Tden);
        fmpz_init(Uden);
        fmpz_init(Vden);

        k = n;
        for (i = 1; (1L << i) < k; i++);
        a = (len_t *) flint_malloc(i * sizeof(len_t));
        a[i = 0] = k;
        while (k >= FLINT_REVERSE_NEWTON_CUTOFF)
            a[++i] = (k = (k + 1) / 2);

        _fmpq_poly_revert_series_lagrange(Qinv, den, Q, Qden, k);
        _fmpz_vec_zero(Qinv + k, n - k);

        for (i--; i >= 0; i--)
        {
            k = a[i];
            _fmpq_poly_compose_series(T, Tden, Q, Qden, k, Qinv, den, k, k);
            _fmpq_poly_derivative(U, Uden, T, Tden, k); fmpz_zero(U + k - 1);
            fmpz_zero(T + 1);
            _fmpq_poly_div_series(V, Vden, T, Tden, U, Uden, k);
            _fmpq_poly_canonicalise(V, Vden, k);
            _fmpq_poly_derivative(T, Tden, Qinv, den, k);
            _fmpq_poly_mullow(U, Uden, V, Vden, k, T, Tden, k, k);
            _fmpq_poly_sub(Qinv, den, Qinv, den, k, U, Uden, k);
        }

        _fmpq_poly_canonicalise(Qinv, den, n);

        flint_free(a);
        _fmpz_vec_clear(T, n);
        _fmpz_vec_clear(U, n);
        _fmpz_vec_clear(V, n);
        fmpz_clear(Tden);
        fmpz_clear(Uden);
        fmpz_clear(Vden);
    }
}


void
fmpq_poly_revert_series_newton(fmpq_poly_t res,
            const fmpq_poly_t poly, len_t n)
{
    fmpz *copy;
    int alloc;

    if (poly->length < 2 || !fmpz_is_zero(poly->coeffs)
                         || fmpz_is_zero(poly->coeffs + 1))
    {
        printf("Exception (fmpq_poly_revert_series_newton). Input must have \n"
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
        len_t i;
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
        _fmpq_poly_revert_series_newton(res->coeffs,
                res->den, copy, poly->den, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_revert_series_newton(t->coeffs,
                t->den, copy, poly->den, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);

    if (alloc)
        flint_free(copy);
}
