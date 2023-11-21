/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dft.h"

static int use_dft(slong len, slong prec)
{
    slong l2 = len;
    while (l2 >= 16) l2 >>= 1;
    if (l2 < 11)
    {
        while (!(len & 1)) len >>= 1;
        while (len % 3 == 0) len /= 3;
        while (len % 5 == 0) len /= 5;
        while (len % 7 == 0) len /= 7;
        return (len == 1);
    }
    return 0;
}

void acb_dft_convol(acb_ptr w, acb_srcptr f, acb_srcptr g, slong len, slong prec)
{
    if (use_dft(len, prec))
        acb_dft_convol_dft(w, f, g, len, prec);
    else
        acb_dft_convol_rad2(w, f, g, len, prec);
}
