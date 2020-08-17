/*
    Copyright (C) 2016 Fredrik Johansson

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

void
_nmod_poly_inv_series_basecase(mp_ptr Qinv, mp_srcptr Q, slong Qlen, slong n, nmod_t mod)
{
    mp_limb_t q;

    Qlen = FLINT_MIN(Qlen, n);

    q = Q[0];
    if (q != 1) q = n_invmod(q, mod.n);
    Qinv[0] = q;

    if (Qlen == 1)
    {
        _nmod_vec_zero(Qinv + 1, n - 1);
    }
    else
    {
        slong i, j, l;
        int nlimbs;
        mp_limb_t s;

        nlimbs = _nmod_vec_dot_bound_limbs(FLINT_MIN(n, Qlen), mod);

        for (i = 1; i < n; i++)
        {
            l = FLINT_MIN(i, Qlen - 1);
            NMOD_VEC_DOT(s, j, l, Q[j + 1], Qinv[i - 1 - j], mod, nlimbs);

            if (q == 1)
                Qinv[i] = n_negmod(s, mod.n);
            else
                Qinv[i] = n_negmod(n_mulmod2_preinv(s, q, mod.n, mod.ninv), mod.n);
        }
    }
}

void
nmod_poly_inv_series_basecase(nmod_poly_t Qinv, const nmod_poly_t Q, slong n)
{
    slong Qlen = Q->length;

    Qlen = FLINT_MIN(Qlen, n);

    if (Qlen == 0)
    {
        flint_printf("Exception (nmod_poly_inv_series_basecase). Division by zero.\n");
        flint_abort();
    }

    if (Qinv != Q)
    {
        nmod_poly_fit_length(Qinv, n);
        _nmod_poly_inv_series_basecase(Qinv->coeffs, Q->coeffs, Qlen, n, Qinv->mod);
    }
    else
    {
        nmod_poly_t t;
        nmod_poly_init2(t, Qinv->mod.n, n);
        _nmod_poly_inv_series_basecase(t->coeffs, Q->coeffs, Qlen, n, Qinv->mod);
        nmod_poly_swap(Qinv, t);
        nmod_poly_clear(t);
    }

    Qinv->length = n;
    _nmod_poly_normalise(Qinv);
}

