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
qqbar_set(qqbar_t res, const qqbar_t x)
{
    fmpz_poly_set(QQBAR_POLY(res), QQBAR_POLY(x));
    acb_set(QQBAR_ENCLOSURE(res), QQBAR_ENCLOSURE(x));
}

