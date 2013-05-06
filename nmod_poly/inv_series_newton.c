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

    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

#define NMOD_POLY_INV_NEWTON_CUTOFF 400

void 
_nmod_poly_inv_series_newton(mp_ptr Qinv, mp_srcptr Q, long n, nmod_t mod)
{
    if (n < NMOD_POLY_INV_NEWTON_CUTOFF)
    {
        _nmod_poly_inv_series_basecase(Qinv, Q, n, mod);
    }
    else
    {
        long *a, i, m;
        mp_ptr W;

        for (i = 1; (1L << i) < n; i++) ;

        W = flint_malloc(n * sizeof(mp_limb_t) + i * sizeof(long));
        a = (long *) (W + n);

        a[i = 0] = n;
        while (n >= NMOD_POLY_INV_NEWTON_CUTOFF)
            a[++i] = (n = (n + 1) / 2);

        _nmod_poly_inv_series_basecase(Qinv, Q, n, mod);

        for (i--; i >= 0; i--)
        {
            m = n;
            n = a[i];

            _nmod_poly_mullow(W, Q, n, Qinv, m, n, mod);
            _nmod_poly_mullow(Qinv + m, Qinv, m, W + m, n - m, n - m, mod);
            _nmod_vec_neg(Qinv + m, Qinv + m, n - m, mod);
        }

        flint_free(W);
    }
}

void
nmod_poly_inv_series_newton(nmod_poly_t Qinv, const nmod_poly_t Q, long n)
{
    const long Qlen = Q->length;

    mp_ptr q, qinv;

    if (n == 0 || Q->length == 0 || Q->coeffs[0] == 0)
    {
        printf("Exception (nmod_poly_inv_series_newton). Division by zero.\n");
        abort();
    }

    if (Qlen < n)
    {
        q = _nmod_vec_init(n);

        flint_mpn_copyi(q, Q->coeffs, Qlen);
        flint_mpn_zero(q + Qlen, n - Qlen);
    }
    else
    {
        q = Q->coeffs;
    }

    if (Q == Qinv && Qlen >= n)
    {
        qinv = _nmod_vec_init(n);
    }
    else
    {
        nmod_poly_fit_length(Qinv, n);
        qinv = Qinv->coeffs;
    }

    _nmod_poly_inv_series_newton(qinv, q, n, Q->mod);

    if (Q == Qinv && Qlen >= n)
    {
        flint_free(Qinv->coeffs);
        Qinv->coeffs = qinv;
        Qinv->alloc  = n;
        Qinv->length = n;
    }
    else
    {
        Qinv->length = n;
    }
    
    if (Qlen < n)
    {
        _nmod_vec_clear(q);
    }

    _nmod_poly_normalise(Qinv);
}

