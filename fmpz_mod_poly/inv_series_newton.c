/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod_poly.h"

void 
_fmpz_mod_poly_inv_series_newton(fmpz * Qinv, const fmpz * Q, slong n, 
                                 const fmpz_t cinv, const fmpz_t p)
{
    if (n == 1)  /* {Q,1} * cinv == 1 mod (x) */
    {
        fmpz_set(Qinv, cinv);
    }
    else
    {
        const slong alloc = FLINT_MAX(n, 3 * FMPZ_MOD_POLY_INV_NEWTON_CUTOFF);
        slong *a, i, m;
        fmpz *W;

        W = _fmpz_vec_init(alloc);

        for (i = 1; (WORD(1) << i) < n; i++) ;

        a = (slong *) flint_malloc(i * sizeof(slong));
        a[i = 0] = n;
        while (n >= FMPZ_MOD_POLY_INV_NEWTON_CUTOFF)
            a[++i] = (n = (n + 1) / 2);

        /* Base case */
        {
            fmpz *Qrev = W + 2 * FMPZ_MOD_POLY_INV_NEWTON_CUTOFF;

            _fmpz_poly_reverse(Qrev, Q, n, n);
            _fmpz_vec_zero(W, 2*n - 2);
            fmpz_one(W + (2*n - 2));
            _fmpz_mod_poly_div_basecase(Qinv, W, W, 2*n - 1, Qrev, n, cinv, p);
            _fmpz_poly_reverse(Qinv, Qinv, n, n);
        }
        
        for (i--; i >= 0; i--)
        {
            m = n;
            n = a[i];

            _fmpz_mod_poly_mullow(W, Q, n, Qinv, m, p, n);
            _fmpz_mod_poly_mullow(Qinv + m, Qinv, m, W + m, n - m, p, n - m);
            _fmpz_mod_poly_neg(Qinv + m, Qinv + m, n - m, p);
        }

        _fmpz_vec_clear(W, alloc);
        flint_free(a);
    }
}

void fmpz_mod_poly_inv_series_newton(fmpz_mod_poly_t Qinv, 
                    const fmpz_mod_poly_t Q, slong n, const fmpz_mod_ctx_t ctx)
{
    const fmpz *p = fmpz_mod_ctx_modulus(ctx);
    fmpz_t cinv;
    fmpz *Qcopy;
    int Qalloc;

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
        flint_mpn_zero((mp_ptr) Qcopy + i, n - i);
        Qalloc = 1;
    }

    fmpz_init(cinv);
    fmpz_invmod(cinv, Q->coeffs, p);

    if (Qinv != Q)
    {
        fmpz_mod_poly_fit_length(Qinv, n, ctx);
        _fmpz_mod_poly_inv_series_newton(Qinv->coeffs, Qcopy, n, cinv, p);
    }
    else
    {
        fmpz *t = _fmpz_vec_init(n);

        _fmpz_mod_poly_inv_series_newton(t, Qcopy, n, cinv, p);

        _fmpz_vec_clear(Qinv->coeffs, Qinv->alloc);
        Qinv->coeffs = t;
        Qinv->alloc  = n;
        Qinv->length = n;
    }
    _fmpz_mod_poly_set_length(Qinv, n);
    _fmpz_mod_poly_normalise(Qinv);

    if (Qalloc)
        flint_free(Qcopy);
    fmpz_clear(cinv);
}

