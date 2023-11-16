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
dlog_crt_clear(dlog_crt_t t)
{
    int k;
    flint_free(t->expo);
    flint_free(t->crt_coeffs);
    for (k = 0; k < t->num; k++)
        dlog_precomp_clear(t->pre + k);
    flint_free(t->pre);
}
