/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "arb_fmpz_poly.h"

void
_arb_fmpz_poly_evaluate_arb(arb_t res, const fmpz * f, slong len,
                           const arb_t x, slong prec)
{
    if (len >= 6 && len >= 5 + 2500 / (FLINT_MAX(prec, 64) + 64))
    {
        /* todo: improve this tuning? */
        if (prec > 1024)
        {
            slong fbits;
            fbits = _fmpz_vec_max_bits(f, len);
            fbits = FLINT_ABS(fbits);

            if (fbits > prec / 2)
            {
                _arb_fmpz_poly_evaluate_arb_horner(res, f, len, x, prec);
                return;
            }
        }

        _arb_fmpz_poly_evaluate_arb_rectangular(res, f, len, x, prec);
        return;
    }

    _arb_fmpz_poly_evaluate_arb_horner(res, f, len, x, prec);
}

void
arb_fmpz_poly_evaluate_arb(arb_t res, const fmpz_poly_t f, const arb_t a, slong prec)
{
    _arb_fmpz_poly_evaluate_arb(res, f->coeffs, f->length, a, prec);
}

