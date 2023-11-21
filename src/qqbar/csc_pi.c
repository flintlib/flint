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
qqbar_csc_pi(qqbar_t res, slong p, ulong q)
{
    qqbar_sin_pi(res, p, q);
    if (qqbar_is_zero(res))
        return 0;
    qqbar_inv(res, res);
    return 1;
}

