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
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

#define FMPQ_POLY_INV_NEWTON_CUTOFF  32

void 
_fmpq_poly_inv_newton(fmpz * rpoly, fmpz_t rden, 
                      const fmpz * poly, const fmpz_t den, long n)
{
    if (n == 1)
    {
        if (fmpz_sgn(poly) > 0)
        {
            fmpz_set(rpoly, den);
            fmpz_set(rden, poly);
        }
        else
        {
            fmpz_neg(rpoly, den);
            fmpz_neg(rden, poly);
        }
    }
    else
    {
        const long alloc = FLINT_MAX(2 * n, 3 * FMPQ_POLY_INV_NEWTON_CUTOFF);
        long *a, i;
        fmpz *W0, *W1, *W0den, *W1den;

        W0 = _fmpz_vec_init(alloc + 2);
        W1 = W0 + n;
        W0den = W0 + alloc;
        W1den = W0 + alloc + 1;

        for (i = 1; (1L << i) < n; i++) ;

        a = (long *) malloc(i * sizeof(long));
        a[i = 0] = n;
        while (n >= FMPQ_POLY_INV_NEWTON_CUTOFF)
            a[++i] = (n = (n + 1) / 2);

        /* Base case */
        {
            fmpz *rev = W0 + 2 * FMPQ_POLY_INV_NEWTON_CUTOFF;

            _fmpz_poly_reverse(rev, poly, n, n);
            _fmpz_vec_zero(W0, 2*n - 2);
            fmpz_set_ui(W0 + (2*n - 2), 1);
            fmpz_set_ui(W0den, 1);

            _fmpq_poly_div(rpoly, rden, W0, W0den, 2*n - 1, 
                                        rev, den, n);

            _fmpz_poly_reverse(rpoly, rpoly, n, n);
        }

        for (i--; i >= 0; i--)
        {
            long len;

            n = a[i];

            _fmpz_poly_mullow_n(W0, poly, rpoly, n);
            fmpz_mul(W0den, den, rden);
            for (len = n - 1; len >= 0 && !W0[len]; len--) ;
            len++;
            _fmpq_poly_canonicalise(W0, W0den, len);

            fmpz_sub(W0, W0, W0den);

            _fmpz_poly_mullow_n(W1, W0, rpoly, n);
            fmpz_mul(W1den, W0den, rden);
            for (len = n - 1; len >= 0 && !W1[len]; len--) ;
            len++;
            _fmpq_poly_canonicalise(W1, W1den, len);

            _fmpq_poly_sub(rpoly, rden, rpoly, rden, n, W1, W1den, n);
        }

        _fmpz_vec_clear(W0, alloc + 2);
        free(a);
    }
}

void fmpq_poly_inv_newton(fmpq_poly_t res, const fmpq_poly_t poly, long n)
{
    fmpz *copy;
    int alloc;

    if (poly->length >= n)
    {
        copy = poly->coeffs;
        alloc = 0;
    }
    else
    {
        long i;
        copy = (fmpz *) malloc(n * sizeof(fmpz));
        for (i = 0; i < poly->length; i++)
            copy[i] = poly->coeffs[i];
        for ( ; i < n; i++)
            copy[i] = 0;
        alloc = 1;
    }

    if (res != poly)
    {
        fmpq_poly_fit_length(res, n);
        _fmpq_poly_inv_newton(res->coeffs, res->den, copy, poly->den, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_inv_newton(t->coeffs, t->den, copy, poly->den, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(res, n);
    fmpq_poly_canonicalise(res);

    if (alloc)
        free(copy);
}

