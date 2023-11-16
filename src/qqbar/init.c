/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "qqbar.h"

void
qqbar_init(qqbar_t res)
{
    fmpz_poly_init(QQBAR_POLY(res));
    fmpz_poly_set_coeff_si(QQBAR_POLY(res), 1, 1);
    acb_init(QQBAR_ENCLOSURE(res));
}

