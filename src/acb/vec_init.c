/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

acb_ptr
_acb_vec_init(slong n)
{
    slong i;
    acb_ptr v = (acb_ptr) flint_malloc(sizeof(acb_struct) * n);

    for (i = 0; i < n; i++)
        acb_init(v + i);

    return v;
}
