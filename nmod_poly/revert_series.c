/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "flint-impl.h"

void
_nmod_poly_revert_series(ulong_ptr Qinv, ulong_srcptr Q, slong n, nmod_t mod)
{
    _nmod_poly_revert_series_lagrange_fast(Qinv, Q, n, mod);
}

void
nmod_poly_revert_series(nmod_poly_t Qinv, 
                                 const nmod_poly_t Q, slong n)
{
    ulong_ptr Qinv_coeffs, Q_coeffs;
    nmod_poly_t t1;
    slong Qlen;
    
    Qlen = Q->length;

    if (Qlen < 2 || Q->coeffs[0] != 0 || Q->coeffs[1] == 0)
        flint_throw(FLINT_ERROR, "Input must have zero constant and an invertible coefficient of x^1 in nmod_poly_revert_series\n");

    if (Qlen < n)
    {
        Q_coeffs = _nmod_vec_init(n);
        FLINT_MPN_COPYI(Q_coeffs, Q->coeffs, Qlen);
        FLINT_MPN_ZERO(Q_coeffs + Qlen, n - Qlen);
    }
    else
        Q_coeffs = Q->coeffs;

    if (Q == Qinv && Qlen >= n)
    {
        nmod_poly_init2(t1, Q->mod.n, n);
        Qinv_coeffs = t1->coeffs;
    }
    else
    {
        nmod_poly_fit_length(Qinv, n);
        Qinv_coeffs = Qinv->coeffs;
    }

    _nmod_poly_revert_series(Qinv_coeffs, Q_coeffs, n, Q->mod);

    if (Q == Qinv && Qlen >= n)
    {
        nmod_poly_swap(Qinv, t1);
        nmod_poly_clear(t1);
    }
    
    Qinv->length = n;

    if (Qlen < n)
        _nmod_vec_clear(Q_coeffs);

    _nmod_poly_normalise(Qinv);
}
