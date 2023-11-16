/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dft.h"

void
acb_dft_convol_naive(acb_ptr w, acb_srcptr f, acb_srcptr g, slong len, slong prec)
{
    slong x, y;
    for (x = 0; x < len; x ++)
    {
        acb_ptr wx;
        acb_srcptr fx;
        wx = w + x;
        fx = f + x;
        acb_zero(wx);
        for (y = 0; y <= x; y++)
            acb_addmul(wx, fx - y, g + y, prec);
        for (; y < len; y++)
            acb_addmul(wx, fx + (len - y), g + y, prec);
    }
}
