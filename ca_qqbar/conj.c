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
ca_qqbar_conj(ca_qqbar_t res, const ca_qqbar_t x)
{
    fmpz_poly_set(CA_QQBAR_POLY(res), CA_QQBAR_POLY(x));
    acb_conj(CA_QQBAR_ENCLOSURE(res), CA_QQBAR_ENCLOSURE(x));
}

