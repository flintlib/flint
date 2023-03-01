/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_revert_series(arb_ptr Qinv,
    arb_srcptr Q, slong Qlen, slong n, slong prec)
{
    _arb_poly_revert_series_lagrange_fast(Qinv, Q, Qlen, n, prec);
}

void
arb_poly_revert_series(arb_poly_t Qinv,
                                    const arb_poly_t Q, slong n, slong prec)
{
    slong Qlen = Q->length;

    if (Qlen < 2 || !arb_is_zero(Q->coeffs)
                 || arb_contains_zero(Q->coeffs + 1))
    {
        flint_printf("Exception (arb_poly_revert_series). Input must \n"
               "have zero constant term and nonzero coefficient of x^1.\n");
        flint_abort();
    }

    if (Qinv != Q)
    {
        arb_poly_fit_length(Qinv, n);
        _arb_poly_revert_series(Qinv->coeffs, Q->coeffs, Qlen, n, prec);
    }
    else
    {
        arb_poly_t t;
        arb_poly_init2(t, n);
        _arb_poly_revert_series(t->coeffs, Q->coeffs, Qlen, n, prec);
        arb_poly_swap(Qinv, t);
        arb_poly_clear(t);
    }

    _arb_poly_set_length(Qinv, n);
    _arb_poly_normalise(Qinv);
}

