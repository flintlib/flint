/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2011, 2016, 2023 Fredrik Johansson
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "gr_poly.h"

void
_nmod_poly_inv_series_basecase_preinv1(mp_ptr Qinv, mp_srcptr Q, slong Qlen, slong n, mp_limb_t q, nmod_t mod)
{
    Qlen = FLINT_MIN(Qlen, n);

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
                Qinv[i] = nmod_neg(s, mod);
            else
                Qinv[i] = nmod_neg(nmod_mul(s, q, mod), mod);
        }
    }
}

void
_nmod_poly_inv_series_basecase(mp_ptr Qinv, mp_srcptr Q, slong Qlen, slong n, nmod_t mod)
{
    mp_limb_t q;

    q = Q[0];
    if (q != 1)
        q = nmod_inv(q, mod);

    _nmod_poly_inv_series_basecase_preinv1(Qinv, Q, Qlen, n, q, mod);
}

void
_nmod_poly_inv_series(mp_ptr Qinv, mp_srcptr Q, slong Qlen, slong n, nmod_t mod)
{
    Qlen = FLINT_MIN(Qlen, n);

    if (Qlen <= 10)
    {
        _nmod_poly_inv_series_basecase(Qinv, Q, Qlen, n, mod);
    }
    else
    {
        gr_ctx_t ctx;
        _gr_ctx_init_nmod(ctx, &mod);
        GR_MUST_SUCCEED(_gr_poly_inv_series(Qinv, Q, Qlen, n, ctx));
    }
}

void
nmod_poly_inv_series(nmod_poly_t Qinv, const nmod_poly_t Q, slong n)
{
    slong Qlen = Q->length;

    Qlen = FLINT_MIN(Qlen, n);

    if (Qlen == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_inv_series). Division by zero.\n");
    }

    if (Qinv != Q)
    {
        nmod_poly_fit_length(Qinv, n);
        _nmod_poly_inv_series(Qinv->coeffs, Q->coeffs, Qlen, n, Qinv->mod);
    }
    else
    {
        nmod_poly_t t;
        nmod_poly_init2(t, Qinv->mod.n, n);
        _nmod_poly_inv_series(t->coeffs, Q->coeffs, Qlen, n, Qinv->mod);
        nmod_poly_swap(Qinv, t);
        nmod_poly_clear(t);
    }

    Qinv->length = n;
    _nmod_poly_normalise(Qinv);
}

void
nmod_poly_inv_series_basecase(nmod_poly_t Qinv, const nmod_poly_t Q, slong n)
{
    slong Qlen = Q->length;

    Qlen = FLINT_MIN(Qlen, n);

    if (Qlen == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_inv_series_basecase). Division by zero.\n");
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

void
_nmod_poly_inv_series_newton(mp_ptr Qinv, mp_srcptr Q, slong Qlen, slong n, nmod_t mod)
{
    _nmod_poly_inv_series(Qinv, Q, Qlen, n, mod);
}

void
nmod_poly_inv_series_newton(nmod_poly_t Qinv, const nmod_poly_t Q, slong n)
{
    nmod_poly_inv_series(Qinv, Q, n);
}
