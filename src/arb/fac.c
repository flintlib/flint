/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_fac_ui(arb_t x, ulong n, slong prec)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_set_ui(t, n);
    fmpz_add_ui(t, t, 1);
    arb_gamma_fmpz(x, t, prec);
    fmpz_clear(t);
}

