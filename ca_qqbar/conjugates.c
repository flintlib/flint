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
ca_qqbar_conjugates(ca_qqbar_ptr res, const ca_qqbar_t x)
{
    if (ca_qqbar_degree(x) == 1)
        ca_qqbar_set(res, x);
    else
        ca_qqbar_roots_fmpz_poly(res, CA_QQBAR_POLY(x), CA_QQBAR_ROOTS_IRREDUCIBLE);
}

