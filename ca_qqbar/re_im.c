/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

void
ca_qqbar_re_im(ca_qqbar_t res1, ca_qqbar_t res2, const ca_qqbar_t x)
{
    if (res1 == x)
    {
        ca_qqbar_im(res2, x);
        ca_qqbar_re(res1, x);
    }
    else
    {
        ca_qqbar_re(res1, x);
        ca_qqbar_im(res2, x);
    }
}
