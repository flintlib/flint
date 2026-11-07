/*
    Copyright (C) 2009, 2010, 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"

ulong
n_div2_preinv(ulong a, ulong n, ulong ninv)
{
    ulong q;

    FLINT_ASSERT(n != 0);

    (void) n_divrem_preinv(&q, a, n, ninv, flint_clz(n));
    return q;
}
