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


#define FLINT_REVERSE_NEWTON_CUTOFF 10

void
_fmpz_poly_revert_series_newton(fmpz * Qinv, const fmpz * Q, long n)
{
    if (n <= 2)
    {
        _fmpz_vec_set(Qinv, Q, n);
        return;
    }
    else
    {
        long *a, i, k;
        fmpz *T, *U, *V;

        T = _fmpz_vec_init(n);
        U = _fmpz_vec_init(n);
        V = _fmpz_vec_init(n);

        k = n;
        for (i = 1; (1L << i) < k; i++);
        a = (long *) flint_malloc(i * sizeof(long));
        a[i = 0] = k;
        while (k >= FLINT_REVERSE_NEWTON_CUTOFF)
            a[++i] = (k = (k + 1) / 2);

        _fmpz_poly_revert_series_lagrange(Qinv, Q, k);
        _fmpz_vec_zero(Qinv + k, n - k);

        for (i--; i >= 0; i--)
        {
            k = a[i];
            _fmpz_poly_compose_series(T, Q, k, Qinv, k, k);
            _fmpz_poly_derivative(U, T, k); fmpz_zero(U + k - 1);
            fmpz_zero(T + 1);
            _fmpz_poly_div_series(V, T, U, k);
            _fmpz_poly_derivative(T, Qinv, k);
            _fmpz_poly_mullow(U, V, k, T, k, k);
            _fmpz_vec_sub(Qinv, Qinv, U, k);
        }

        flint_free(a);
        _fmpz_vec_clear(T, n);
        _fmpz_vec_clear(U, n);
        _fmpz_vec_clear(V, n);
    }
}

void
fmpz_poly_revert_series_newton(fmpz_poly_t Qinv, const fmpz_poly_t Q, long n)
{
    fmpz *Qcopy;
    int Qalloc;
    long Qlen = Q->length;

    if (Qlen < 2 || !fmpz_is_zero(Q->coeffs) || !fmpz_is_pm1(Q->coeffs + 1))
    {
        printf("Exception (fmpz_poly_revert_series_newton). Input must have \n"
               "zero constant term and +1 or -1 as coefficient of x^1.\n");
        abort();
    }

    if (Qlen >= n)
    {
        Qcopy = Q->coeffs;
        Qalloc = 0;
    }
    else
    {
        long i;
        Qcopy = (fmpz *) flint_malloc(n * sizeof(fmpz));
        for (i = 0; i < Qlen; i++)
            Qcopy[i] = Q->coeffs[i];
        for ( ; i < n; i++)
            Qcopy[i] = 0;
        Qalloc = 1;
    }

    if (Qinv != Q)
    {
        fmpz_poly_fit_length(Qinv, n);
        _fmpz_poly_revert_series_newton(Qinv->coeffs, Qcopy, n);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        _fmpz_poly_revert_series_newton(t->coeffs, Qcopy, n);
        fmpz_poly_swap(Qinv, t);
        fmpz_poly_clear(t);
    }

    _fmpz_poly_set_length(Qinv, n);
    _fmpz_poly_normalise(Qinv);

    if (Qalloc)
        flint_free(Qcopy);
}

