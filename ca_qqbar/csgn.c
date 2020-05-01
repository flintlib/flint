/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

int
ca_qqbar_csgn(const ca_qqbar_t x)
{
    int re = ca_qqbar_real_sgn(x);

    if (re != 0)
        return re;

    return ca_qqbar_imag_sgn(x);
}
