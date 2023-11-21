/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qqbar.h"

int
qqbar_csgn(const qqbar_t x)
{
    int re = qqbar_sgn_re(x);

    if (re != 0)
        return re;

    return qqbar_sgn_im(x);
}
