/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

void
mag_log_ui(mag_t t, ulong n)
{
    if (n <= 1)
    {
        if (n == 1)
            mag_zero(t);
        else
            mag_inf(t);
    }
    else
    {
        mag_set_ui(t, n-1);
        mag_log1p(t, t);
    }
}

