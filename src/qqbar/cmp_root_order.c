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
qqbar_cmp_root_order(const qqbar_t x, const qqbar_t y)
{
    int xreal, yreal, cmp;

    xreal = qqbar_is_real(x);
    yreal = qqbar_is_real(y);

    if (xreal != yreal)
        return xreal ? -1 : 1;

    cmp = qqbar_cmp_re(x, y);
    if (cmp != 0)
        return -cmp;

    cmp = qqbar_cmpabs_im(x, y);
    if (cmp != 0)
        return cmp;

    return qqbar_sgn_im(y);
}
