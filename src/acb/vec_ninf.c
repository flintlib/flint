/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
_acb_vec_ninf(arb_t ninf, acb_srcptr vec, slong len, slong prec)
{
    arb_t abs;
    slong i;
    
    arb_init(abs);
    arb_zero(ninf);

    for (i = 0; i < len; i++)
    {
        acb_abs(abs, vec + i, prec);
        arb_max(ninf, ninf, abs, prec);
    }

    arb_clear(abs);
}
