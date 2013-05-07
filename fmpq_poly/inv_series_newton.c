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

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

#define FMPQ_POLY_INV_NEWTON_CUTOFF  32

void 
_fmpq_poly_inv_series_newton(fmpz * Qinv, fmpz_t Qinvden, 
                             const fmpz * Q, const fmpz_t Qden, long n)
{
    if (n == 1)
    {
        if (fmpz_sgn(Q) > 0)
        {
            fmpz_set(Qinv, Qden);
            fmpz_set(Qinvden, Q);
        }
        else
        {
            fmpz_neg(Qinv, Qden);
            fmpz_neg(Qinvden, Q);
        }
    }
    else
    {
        const long alloc = FLINT_MAX(n, 3 * FMPQ_POLY_INV_NEWTON_CUTOFF);
        long *a, i, m;
        fmpz *W, *Wden;

        W = _fmpz_vec_init(alloc + 1);
        Wden = W + alloc;

        for (i = 1; (1L << i) < n; i++) ;

        a = (long *) flint_malloc(i * sizeof(long));
        a[i = 0] = n;
        while (n >= FMPQ_POLY_INV_NEWTON_CUTOFF)
            a[++i] = (n = (n + 1) / 2);

        /* Base case */
        {
            fmpz *rev = W + 2 * FMPQ_POLY_INV_NEWTON_CUTOFF;

            _fmpz_poly_reverse(rev, Q, n, n);
            _fmpz_vec_zero(W, 2*n - 2);
            fmpz_one(W + (2*n - 2));
            fmpz_one(Wden);

            _fmpq_poly_div(Qinv, Qinvden, W, Wden, 2*n - 1, rev, Qden, n);
            _fmpq_poly_canonicalise(Qinv, Qinvden, n);

            _fmpz_poly_reverse(Qinv, Qinv, n, n);
        }

        for (i--; i >= 0; i--)
        {
            m = n;
            n = a[i];

            _fmpz_poly_mullow(W, Q, n, Qinv, m, n);
            fmpz_mul(Wden, Qden, Qinvden);

            _fmpz_poly_mullow(Qinv + m, Qinv, m, W + m, n - m, n - m);
            fmpz_mul(Qinvden, Qinvden, Wden);
            _fmpz_vec_scalar_mul_fmpz(Qinv, Qinv, m, Wden);

            _fmpz_vec_neg(Qinv + m, Qinv + m, n - m);

            _fmpq_poly_canonicalise(Qinv, Qinvden, n);
        }

        _fmpz_vec_clear(W, alloc + 1);
        flint_free(a);
    }
}

void fmpq_poly_inv_series_newton(fmpq_poly_t Qinv, const fmpq_poly_t Q, long n)
{
    fmpz *copy;
    int alloc;

    if (Q->length >= n)
    {
        copy = Q->coeffs;
        alloc = 0;
    }
    else
    {
        long i;
        copy = (fmpz *) flint_malloc(n * sizeof(fmpz));
        for (i = 0; i < Q->length; i++)
            copy[i] = Q->coeffs[i];
        for ( ; i < n; i++)
            copy[i] = 0;
        alloc = 1;
    }

    if (Qinv != Q)
    {
        fmpq_poly_fit_length(Qinv, n);
        _fmpq_poly_inv_series_newton(Qinv->coeffs, Qinv->den, copy, Q->den, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_inv_series_newton(t->coeffs, t->den, copy, Q->den, n);
        fmpq_poly_swap(Qinv, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(Qinv, n);
    fmpq_poly_canonicalise(Qinv);

    if (alloc)
        flint_free(copy);
}

