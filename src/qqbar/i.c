/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

void
qqbar_i(qqbar_t res)
{
    fmpz_poly_zero(QQBAR_POLY(res));
    fmpz_poly_set_coeff_si(QQBAR_POLY(res), 2, 1);
    fmpz_poly_set_coeff_si(QQBAR_POLY(res), 0, 1);
    acb_onei(QQBAR_ENCLOSURE(res));
}

