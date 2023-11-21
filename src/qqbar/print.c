/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

void
qqbar_print(const qqbar_t x)
{
    slong i, d;
    d = qqbar_degree(x);

    flint_printf("deg %wd [", qqbar_degree(x));
    for (i = 0; i <= d; i++)
    {
        fmpz_print(QQBAR_COEFFS(x) + i);
        if (i < d)
            flint_printf(", ");
    }
    flint_printf("] ");
    acb_printn(QQBAR_ENCLOSURE(x), FLINT_MAX(6, FLINT_MIN(acb_rel_accuracy_bits(QQBAR_ENCLOSURE(x)),
        acb_bits(QQBAR_ENCLOSURE(x)))), 0);
}

