/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

arb_ptr
_arb_vec_init(slong n)
{
    slong i;
    arb_ptr v = (arb_ptr) flint_malloc(sizeof(arb_struct) * n);

    for (i = 0; i < n; i++)
        arb_init(v + i);

    return v;
}
