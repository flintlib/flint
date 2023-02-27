/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2011, 2016 Fredrik Johansson
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

#define MULLOW(z, x, xn, y, yn, nn, mod) \
    if ((xn) >= (yn)) \
        _nmod_poly_mullow(z, x, xn, y, yn, nn, mod); \
    else \
        _nmod_poly_mullow(z, y, yn, x, xn, nn, mod); \

void
_nmod_poly_inv_series_newton(mp_ptr Qinv, mp_srcptr Q, slong Qlen, slong n, nmod_t mod)
{
    slong cutoff;

    Qlen = FLINT_MIN(Qlen, n);

    if (Qlen < 16 || mod.n <= 3)
        cutoff = 16;
    else
        cutoff = 25 * FLINT_BIT_COUNT(mod.n);

    if (Qlen < cutoff)
    {
        _nmod_poly_inv_series_basecase(Qinv, Q, Qlen, n, mod);
    }
    else
    {
        slong *a, i, m, Qnlen, Wlen, W2len;
        mp_ptr W;

        for (i = 1; (WORD(1) << i) < n; i++) ;

        W = flint_malloc(n * sizeof(mp_limb_t) + i * sizeof(slong));
        a = (slong *) (W + n);

        a[i = 0] = n;
        while (n >= cutoff)
            a[++i] = (n = (n + 1) / 2);

        _nmod_poly_inv_series_basecase(Qinv, Q, Qlen, n, mod);

        for (i--; i >= 0; i--)
        {
            m = n;
            n = a[i];

            Qnlen = FLINT_MIN(Qlen, n);
            Wlen = FLINT_MIN(Qnlen + m - 1, n);
            W2len = Wlen - m;
            MULLOW(W, Q, Qnlen, Qinv, m, Wlen, mod);
            MULLOW(Qinv + m, Qinv, m, W + m, W2len, n - m, mod);
            _nmod_vec_neg(Qinv + m, Qinv + m, n - m, mod);
        }

        flint_free(W);
    }
}

void
nmod_poly_inv_series_newton(nmod_poly_t Qinv, const nmod_poly_t Q, slong n)
{
    slong Qlen = Q->length;

    Qlen = FLINT_MIN(Qlen, n);

    if (Qlen == 0)
    {
        flint_printf("Exception (nmod_poly_inv_series_newton). Division by zero.\n");
        flint_abort();
    }

    if (Qinv != Q)
    {
        nmod_poly_fit_length(Qinv, n);
        _nmod_poly_inv_series_newton(Qinv->coeffs, Q->coeffs, Qlen, n, Qinv->mod);
    }
    else
    {
        nmod_poly_t t;
        nmod_poly_init2(t, Qinv->mod.n, n);
        _nmod_poly_inv_series_newton(t->coeffs, Q->coeffs, Qlen, n, Qinv->mod);
        nmod_poly_swap(Qinv, t);
        nmod_poly_clear(t);
    }

    Qinv->length = n;
    _nmod_poly_normalise(Qinv);
}

