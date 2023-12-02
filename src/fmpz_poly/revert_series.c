/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011-2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"
#include "fmpz.h"
#include "fmpz_poly.h"

void
_fmpz_poly_revert_series(fmpz * Qinv, const fmpz * Q, slong Qlen, slong n)
{
    gr_ctx_t ctx;
    gr_ctx_init_fmpz(ctx);
    GR_MUST_SUCCEED(_gr_poly_revert_series(Qinv, Q, Qlen, n, ctx));
}

void
fmpz_poly_revert_series(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n)
{
    slong Qlen = Q->length;

    if (Qlen < 2 || !fmpz_is_zero(Q->coeffs) || !fmpz_is_pm1(Q->coeffs + 1))
    {
        flint_throw(FLINT_ERROR, "(fmpz_poly_revert_series): "
                "Input must have zero constant term and +1 or -1 as coefficient of x^1\n.");
    }

    if (Qinv != Q)
    {
        fmpz_poly_fit_length(Qinv, n);
        _fmpz_poly_revert_series(Qinv->coeffs, Q->coeffs, Qlen, n);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        _fmpz_poly_revert_series(t->coeffs, Q->coeffs, Qlen, n);
        fmpz_poly_swap(Qinv, t);
        fmpz_poly_clear(t);
    }

    _fmpz_poly_set_length(Qinv, n);
    _fmpz_poly_normalise(Qinv);
}
