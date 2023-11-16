/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "dlog.h"

void
dlog_vec_fill(ulong *v, ulong nv, ulong x)
{
    ulong k;
    for (k = 0; k < nv; k++)
        v[k] = x;
}
