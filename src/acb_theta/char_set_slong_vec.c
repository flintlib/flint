/*
    Copyright (C) 2025 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

ulong
acb_theta_char_set_slong_vec(const slong * vec, slong len)
{
    slong k;
    ulong ch = 0;

    for (k = 0; k < len; k++)
    {
        ch = ch << 1;
        ch += (vec[k] % 2);
    }
    return ch;
}
