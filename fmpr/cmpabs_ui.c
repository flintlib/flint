/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpr.h"

int
fmpr_cmpabs_ui(const fmpr_t x, ulong y)
{
    fmpr_t t;
    int res;
    fmpr_init(t);
    fmpr_set_ui(t, y);
    res = fmpr_cmpabs(x, t);
    fmpr_clear(t);
    return res;
}

