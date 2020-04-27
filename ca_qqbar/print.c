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
ca_qqbar_print(const ca_qqbar_t x)
{
    flint_printf("ca_qqbar with poly = {");
    fmpz_poly_print(CA_QQBAR_POLY(x));
    flint_printf("} and enclosure = ");
    acb_printn(CA_QQBAR_ENCLOSURE(x), FLINT_MAX(6, FLINT_MIN(acb_rel_accuracy_bits(CA_QQBAR_ENCLOSURE(x)),
        acb_bits(CA_QQBAR_ENCLOSURE(x)))), 0);
}

