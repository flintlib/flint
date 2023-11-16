/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_evaluate2(acb_t y, acb_t z, acb_srcptr f, slong len, const acb_t x, slong prec)
{
    if ((prec >= 1024) && (len >= 5 + 20000 / prec))
    {
        slong fbits;

        fbits = _acb_vec_bits(f, len);

        if (fbits <= prec / 2)
        {
            _acb_poly_evaluate2_rectangular(y, z, f, len, x, prec);
            return;
        }
    }

    _acb_poly_evaluate2_horner(y, z, f, len, x, prec);
}

void
acb_poly_evaluate2(acb_t r, acb_t s, const acb_poly_t f, const acb_t a, slong prec)
{
    _acb_poly_evaluate2(r, s, f->coeffs, f->length, a, prec);
}

