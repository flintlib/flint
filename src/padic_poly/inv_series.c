/*
    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2023

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "padic.h"
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
        const slong alloc = FLINT_MAX(n / 2, 3);
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

            _fmpz_poly_mulmid(W, Q, n, Qinv, m, m, n);
            _fmpz_vec_scalar_mod_fmpz(W, W, n - m, p);
            _fmpz_poly_mullow(Qinv + m, Qinv, m, W, n - m, n - m);
            _fmpz_vec_neg(Qinv + m, Qinv + m, n - m);
            _fmpz_vec_scalar_mod_fmpz(Qinv + m, Qinv + m, n - m, p);
        }

        _fmpz_vec_clear(W, alloc);
    }
}

void padic_poly_inv_series(padic_poly_t g, const padic_poly_t f, slong n,
                           const padic_ctx_t ctx)
{
    fmpz_t cinv;

    fmpz_t pow;
    int palloc;

    fmpz *Qcopy;
    int Qalloc;

    if (f->length == 0 || fmpz_is_zero(f->coeffs + 0))
    {
        flint_throw(FLINT_ERROR, "Exception (padic_poly_inv_series):  Constant term is zero.\n");
    }
    if (fmpz_divisible(f->coeffs + 0, ctx->p))
    {
        flint_throw(FLINT_ERROR, "Exception (padic_poly_inv_series):\nValuation of constant term is not minimal.\n");
    }

    if (- f->val >= g->N)
    {
        padic_poly_zero(g);
        return;
    }

    if (f->length >= n)
    {
        Qcopy = f->coeffs;
        Qalloc = 0;
    }
    else
    {
        slong i;

        Qcopy = (fmpz *) flint_malloc(n * sizeof(fmpz));
        for (i = 0; i < f->length; i++)
            Qcopy[i] = f->coeffs[i];
        mpn_zero((nn_ptr) Qcopy + i, n - i);
        Qalloc = 1;
    }

    fmpz_init(cinv);
    fmpz_init(pow);

    _padic_inv(cinv, f->coeffs, ctx->p, g->N + f->val);
    palloc = _padic_ctx_pow_ui(pow, g->N + f->val, ctx);

    if (g != f)
    {
        padic_poly_fit_length(g, n);

        /* fails: _fmpz_vec_scalar_mod_fmpz(Qcopy, Qcopy, n, pow); */
        _fmpz_mod_poly_inv_series_pure_newton(g->coeffs, Qcopy, n, cinv, pow);
    }
    else
    {
        fmpz *t = _fmpz_vec_init(n);

        /* fails: _fmpz_vec_scalar_mod_fmpz(Qcopy, Qcopy, n, pow); */
        _fmpz_mod_poly_inv_series_pure_newton(t, Qcopy, n, cinv, pow);
        _fmpz_vec_clear(g->coeffs, g->alloc);
        g->coeffs = t;
        g->alloc  = n;
        g->length = n;
    }

    g->val = - f->val;

    _padic_poly_set_length(g, n);
    _padic_poly_normalise(g);

    fmpz_clear(cinv);
    if (palloc)
        fmpz_clear(pow);
    if (Qalloc)
        flint_free(Qcopy);
}
