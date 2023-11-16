/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

void
arb_hypgeom_rising_ui_jet(arb_ptr res, const arb_t x, ulong n, slong len, slong prec)
{
    if (len == 1)
    {
        arb_hypgeom_rising_ui_rec(res, x, n, prec);
    }
    else if (n <= 7)
    {
        arb_hypgeom_rising_ui_jet_powsum(res, x, n, len, prec);
    }
    else if (len == 2)
    {
        if (n <= 30 || arb_bits(x) >= prec / 128)
            arb_hypgeom_rising_ui_jet_rs(res, x, n, 0, len, prec);
        else
            arb_hypgeom_rising_ui_jet_bs(res, x, n, len, prec);
    }
    else
    {
        if (n <= 20 || (n <= 200 && prec > 400 * n && arb_bits(x) >= prec / 4))
        {
            arb_hypgeom_rising_ui_jet_powsum(res, x, n, len, prec);
        }
        else if (len >= 64 || (arb_bits(x) + 1 < prec / 1024 && n >= 32))
        {
            arb_hypgeom_rising_ui_jet_bs(res, x, n, len, prec);
        }
        else
        {
            arb_hypgeom_rising_ui_jet_rs(res, x, n, 0, len, prec);
        }
    }
}

