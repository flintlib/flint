/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_pow_fmpz(arb_t y, const arb_t b, const fmpz_t e, slong prec)
{
    arb_pow_fmpz_binexp(y, b, e, prec);
}

void
arb_pow_ui(arb_t y, const arb_t b, ulong e, slong prec)
{
    fmpz_t f;
    fmpz_init_set_ui(f, e);
    arb_pow_fmpz(y, b, f, prec);
    fmpz_clear(f);
}
