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
ca_qqbar_init(ca_qqbar_t res)
{
    fmpz_poly_init(CA_QQBAR_POLY(res));
    fmpz_poly_set_coeff_si(CA_QQBAR_POLY(res), 1, 1);
    acb_init(CA_QQBAR_ENCLOSURE(res));
}

