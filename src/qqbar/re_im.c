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
qqbar_re_im(qqbar_t res1, qqbar_t res2, const qqbar_t x)
{
    if (res1 == x)
    {
        qqbar_im(res2, x);
        qqbar_re(res1, x);
    }
    else
    {
        qqbar_re(res1, x);
        qqbar_im(res2, x);
    }
}
