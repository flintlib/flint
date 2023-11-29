/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"
#include "acb_poly.h"

void
_acb_poly_revert_series(acb_ptr Qinv,
    acb_srcptr Q, slong Qlen, slong n, slong prec)
{
    gr_ctx_t ctx;
    gr_ctx_init_complex_acb(ctx, prec);
    if (_gr_poly_revert_series(Qinv, Q, Qlen, n, ctx) != GR_SUCCESS)
        _acb_vec_indeterminate(Qinv, n);
}

void
acb_poly_revert_series(acb_poly_t Qinv,
                                    const acb_poly_t Q, slong n, slong prec)
{
    slong Qlen = Q->length;

    if (Qlen < 2 || !acb_is_zero(Q->coeffs)
                 || acb_contains_zero(Q->coeffs + 1))
    {
        flint_throw(FLINT_ERROR, "(acb_poly_revert_series): Input must \n"
               "have zero constant term and nonzero coefficient of x^1.\n");
    }

    if (Qinv != Q)
    {
        acb_poly_fit_length(Qinv, n);
        _acb_poly_revert_series(Qinv->coeffs, Q->coeffs, Qlen, n, prec);
    }
    else
    {
        acb_poly_t t;
        acb_poly_init2(t, n);
        _acb_poly_revert_series(t->coeffs, Q->coeffs, Qlen, n, prec);
        acb_poly_swap(Qinv, t);
        acb_poly_clear(t);
    }

    _acb_poly_set_length(Qinv, n);
    _acb_poly_normalise(Qinv);
}
