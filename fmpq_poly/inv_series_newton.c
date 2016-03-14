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
    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

#define FMPQ_POLY_INV_NEWTON_CUTOFF 24

/* Requires 2*min(Qlen,n) + n - 1 < 3n coefficients of scratch space in W */
static void
_fmpq_poly_inv_series_basecase_rev(fmpz * Qinv, fmpz_t Qinvden,
    fmpz * W, fmpz_t Wden,
    const fmpz * Q, const fmpz_t Qden, slong Qlen, slong n)
{
    slong Wlen;
    fmpz *Qrev;

    Qlen = FLINT_MIN(Qlen, n);
    Wlen = n + Qlen - 1;
    Qrev = W + Wlen;

    _fmpz_poly_reverse(Qrev, Q, Qlen, Qlen);
    _fmpz_vec_zero(W, Wlen - 1);
    fmpz_one(W + Wlen - 1);
    fmpz_one(Wden);

    _fmpq_poly_div(Qinv, Qinvden, W, Wden, Wlen, Qrev, Qden, Qlen, NULL);

    _fmpq_poly_canonicalise(Qinv, Qinvden, n);
    _fmpz_poly_reverse(Qinv, Qinv, n, n);
}


#define MULLOW(z, x, xn, y, yn, nn) \
    if ((xn) >= (yn)) \
        _fmpz_poly_mullow(z, x, xn, y, yn, nn); \
    else \
        _fmpz_poly_mullow(z, y, yn, x, xn, nn); \

void 
_fmpq_poly_inv_series_newton(fmpz * Qinv, fmpz_t Qinvden, 
                     const fmpz * Q, const fmpz_t Qden, slong Qlen, slong n)
{
    if (fmpz_is_pm1(Q) && fmpz_is_one(Qden))
    {
        _fmpz_poly_inv_series(Qinv, Q, Qlen, n);
        fmpz_one(Qinvden);
        return;
    }

    Qlen = FLINT_MIN(Qlen, n);

    if (Qlen == 1)
    {
        fmpz_set(Qinv, Qden);
        fmpz_set(Qinvden, Q);
        _fmpq_canonicalise(Qinv, Qinvden);
        _fmpz_vec_zero(Qinv + 1, n - 1);
    }
    else
    {
        slong alloc, Qnlen, Wlen, W2len;
        fmpz * W;
        fmpz_t Wden;

        alloc = FLINT_MAX(n, 3 * FMPQ_POLY_INV_NEWTON_CUTOFF);
        W = _fmpz_vec_init(alloc);
        fmpz_init(Wden);

        FLINT_NEWTON_INIT(FMPQ_POLY_INV_NEWTON_CUTOFF, n)

        FLINT_NEWTON_BASECASE(n)
        _fmpq_poly_inv_series_basecase_rev(Qinv, Qinvden, W, Wden, Q, Qden, Qlen, n);
        FLINT_NEWTON_END_BASECASE

        FLINT_NEWTON_LOOP(m, n)

        Qnlen = FLINT_MIN(Qlen, n);
        Wlen = FLINT_MIN(Qnlen + m - 1, n);
        W2len = Wlen - m;

        MULLOW(W, Q, Qnlen, Qinv, m, Wlen);
        fmpz_mul(Wden, Qden, Qinvden);

        MULLOW(Qinv + m, Qinv, m, W + m, W2len, n - m);
        fmpz_mul(Qinvden, Qinvden, Wden);

        _fmpz_vec_scalar_mul_fmpz(Qinv, Qinv, m, Wden);
        _fmpz_vec_neg(Qinv + m, Qinv + m, n - m);
        _fmpq_poly_canonicalise(Qinv, Qinvden, n);

        FLINT_NEWTON_END_LOOP

        FLINT_NEWTON_END

        _fmpz_vec_clear(W, alloc);
        fmpz_clear(Wden);
    }
}

void
fmpq_poly_inv_series_newton(fmpq_poly_t Qinv, const fmpq_poly_t Q, slong n)
{
    slong Qlen = Q->length;

    Qlen = FLINT_MIN(Qlen, n);

    if (Qlen == 0)
    {
        flint_printf("Exception (fmpq_poly_inv_series_newton). Division by zero.\n");
        flint_abort();
    }

    if (Qinv != Q)
    {
        fmpq_poly_fit_length(Qinv, n);
        _fmpq_poly_inv_series_newton(Qinv->coeffs, Qinv->den, Q->coeffs, Q->den, Qlen, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_inv_series_newton(t->coeffs, t->den, Q->coeffs, Q->den, Qlen, n);
        fmpq_poly_swap(Qinv, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(Qinv, n);
    _fmpq_poly_normalise(Qinv);
}

