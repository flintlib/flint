/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_fmpz_poly.h"

void
_arb_fmpz_poly_evaluate_acb(acb_t res, const fmpz * f, slong len,
                           const acb_t x, slong prec)
{
    if (acb_is_real(x))
    {
        _arb_fmpz_poly_evaluate_arb(acb_realref(res), f, len, acb_realref(x), prec);
        arb_zero(acb_imagref(res));
    }
    else
    {
        _arb_fmpz_poly_evaluate_acb_rectangular(res, f, len, x, prec);
    }
}

void
arb_fmpz_poly_evaluate_acb(acb_t res, const fmpz_poly_t f, const acb_t a, slong prec)
{
    _arb_fmpz_poly_evaluate_acb(res, f->coeffs, f->length, a, prec);
}

