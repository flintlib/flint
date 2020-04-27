/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

void
ca_qqbar_swap(ca_qqbar_t x, ca_qqbar_t y)
{
    fmpz_poly_swap(CA_QQBAR_POLY(x), CA_QQBAR_POLY(y));
    acb_swap(CA_QQBAR_ENCLOSURE(x), CA_QQBAR_ENCLOSURE(y));
}

