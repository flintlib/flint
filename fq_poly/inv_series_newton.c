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

    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

******************************************************************************/

#include "fq_poly.h"

#define FQ_POLY_INV_NEWTON_CUTOFF  64

void 
_fq_poly_inv_series_newton(fq_struct * Qinv, const fq_struct * Q, slong n, 
                           const fq_t cinv, const fq_ctx_t ctx)
{
    if (n == 1)  /* {Q,1} x* cinv == 1 mod (x) */
    {
        fq_set(Qinv, cinv, ctx);
    }
    else
    {
        const slong alloc = FLINT_MAX(n, 3 * FQ_POLY_INV_NEWTON_CUTOFF);
        slong *a, i, m;
        fq_struct *W;

        W = _fq_vec_init(alloc, ctx);

        for (i = 1; (WORD(1) << i) < n; i++) ;

        a = (slong *) flint_malloc(i * sizeof(slong));
        a[i = 0] = n;
        while (n >= FQ_POLY_INV_NEWTON_CUTOFF)
            a[++i] = (n = (n + 1) / 2);

        /* Base case */
        {
            fq_struct *Qrev = W + 2 * FQ_POLY_INV_NEWTON_CUTOFF;

            _fq_poly_reverse(Qrev, Q, n, n, ctx);
            _fq_vec_zero(W, 2*n - 2, ctx);
            fq_one(W + (2*n - 2), ctx);
            _fq_poly_div_basecase(Qinv, W, W, 2*n - 1, Qrev, n, cinv, ctx);
            _fq_poly_reverse(Qinv, Qinv, n, n, ctx);
        }
        
        for (i--; i >= 0; i--)
        {
            m = n;
            n = a[i];

            _fq_poly_mullow(W, Q, n, Qinv, m, n, ctx);
            _fq_poly_mullow(Qinv + m, Qinv, m, W + m, n - m, n - m, ctx);
            _fq_poly_neg(Qinv + m, Qinv + m, n - m, ctx);
        }

        _fq_vec_clear(W, alloc, ctx);
        flint_free(a);
    }
}

void fq_poly_inv_series_newton(fq_poly_t Qinv, const fq_poly_t Q, slong n,
                               const fq_ctx_t ctx)
{
    fq_t cinv;
    fq_struct *Qcopy;
    int Qalloc;

    if (Q->length >= n)
    {
        Qcopy = Q->coeffs;
        Qalloc = 0;
    }
    else
    {
        Qcopy = _fq_vec_init(n, ctx);
        _fq_vec_set(Qcopy, Q->coeffs, Q->length, ctx);
        Qalloc = 1;
    }

    fq_init(cinv, ctx);
    fq_inv(cinv, Q->coeffs, ctx);

    if (Qinv != Q)
    {
        fq_poly_fit_length(Qinv, n, ctx);
        _fq_poly_inv_series_newton(Qinv->coeffs, Qcopy, n, cinv, ctx);
    }
    else
    {
        fq_struct *t = _fq_vec_init(n, ctx);

        _fq_poly_inv_series_newton(t, Qcopy, n, cinv, ctx);

        _fq_vec_clear(Qinv->coeffs, Qinv->alloc, ctx);
        Qinv->coeffs = t;
        Qinv->alloc  = n;
        Qinv->length = n;
    }
    _fq_poly_set_length(Qinv, n, ctx);
    _fq_poly_normalise(Qinv, ctx);

    if (Qalloc)
        _fq_vec_clear(Qcopy, n, ctx);
    fq_clear(cinv, ctx);
}

