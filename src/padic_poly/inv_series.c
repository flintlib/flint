/*
    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2023

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "padic_poly.h"

/*  The tests fail if one reduces the coefficients modulo N before
    doing the calculation in Z/NZ. Therefore _fmpz_mod_poly_inv_series
    cannot be used. To do: explain why this particular algorithm,
    given non-normalized input, gives the expected result. */

static void
_fmpz_mod_poly_inv_series_pure_newton(fmpz * Qinv, const fmpz * Q, slong n,
                                 const fmpz_t cinv, const fmpz_t p)
{
    if (n == 1)
    {
        fmpz_set(Qinv, cinv);
    }
    else
    {
        const slong alloc = FLINT_MAX(n, 3);
        slong a[FLINT_BITS], i, m;
        fmpz *W;

        W = _fmpz_vec_init(alloc);

        for (i = 1; (WORD(1) << i) < n; i++) ;

        a[i = 0] = n;
        while (n >= 2)
            a[++i] = (n = (n + 1) / 2);

        /* Base case */
        fmpz_set(Qinv, cinv);

        for (i--; i >= 0; i--)
        {
            m = n;
            n = a[i];

            _fmpz_poly_mullow(W, Q, n, Qinv, m, n);
            _fmpz_vec_scalar_mod_fmpz(W, W, n, p);
            _fmpz_poly_mullow(Qinv + m, Qinv, m, W + m, n - m, n - m);
            _fmpz_vec_neg(Qinv + m, Qinv + m, n - m);
            _fmpz_vec_scalar_mod_fmpz(Qinv + m, Qinv + m, n - m, p);
        }

        _fmpz_vec_clear(W, alloc);
    }
}

void padic_poly_inv_series(padic_poly_t Qinv, const padic_poly_t Q, slong n,
                           const padic_ctx_t ctx)
{
    fmpz_t cinv;

    fmpz_t pow;
    int palloc;

    fmpz *Qcopy;
    int Qalloc;

    if (Q->length == 0 || fmpz_is_zero(Q->coeffs + 0))
    {
        flint_throw(FLINT_ERROR, "Exception (padic_poly_inv_series):  Constant term is zero.\n");
    }
    if (fmpz_divisible(Q->coeffs + 0, ctx->p))
    {
        flint_throw(FLINT_ERROR, "Exception (padic_poly_inv_series):\nValuation of constant term is not minimal.\n");
    }

    if (- Q->val >= Qinv->N)
    {
        padic_poly_zero(Qinv);
        return;
    }

    if (Q->length >= n)
    {
        Qcopy = Q->coeffs;
        Qalloc = 0;
    }
    else
    {
        slong i;

        Qcopy = (fmpz *) flint_malloc(n * sizeof(fmpz));
        for (i = 0; i < Q->length; i++)
            Qcopy[i] = Q->coeffs[i];
        mpn_zero((mp_ptr) Qcopy + i, n - i);
        Qalloc = 1;
    }

    fmpz_init(cinv);
    fmpz_init(pow);

    _padic_inv(cinv, Q->coeffs, ctx->p, Qinv->N + Q->val);
    palloc = _padic_ctx_pow_ui(pow, Qinv->N + Q->val, ctx);

    if (Qinv != Q)
    {
        padic_poly_fit_length(Qinv, n);

        /* fails: _fmpz_vec_scalar_mod_fmpz(Qcopy, Qcopy, n, pow); */
        _fmpz_mod_poly_inv_series_pure_newton(Qinv->coeffs, Qcopy, n, cinv, pow);
    }
    else
    {
        fmpz *t = _fmpz_vec_init(n);

        /* fails: _fmpz_vec_scalar_mod_fmpz(Qcopy, Qcopy, n, pow); */
        _fmpz_mod_poly_inv_series_pure_newton(t, Qcopy, n, cinv, pow);
        _fmpz_vec_clear(Qinv->coeffs, Qinv->alloc);
        Qinv->coeffs = t;
        Qinv->alloc  = n;
        Qinv->length = n;
    }

    Qinv->val = - Q->val;

    _padic_poly_set_length(Qinv, n);
    _padic_poly_normalise(Qinv);

    fmpz_clear(cinv);
    if (palloc)
        fmpz_clear(pow);
    if (Qalloc)
        flint_free(Qcopy);
}

