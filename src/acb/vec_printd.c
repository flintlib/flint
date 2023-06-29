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
_acb_vec_printd(acb_srcptr vec, slong len, slong digits)
{
    slong i;
    flint_printf("[");
    for (i = 0; i < len; i++)
    {
        if (i > 0)
        {
            flint_printf(" ");
        }
        acb_printd(vec + i, digits);
        if (i < len - 1)
        {
            flint_printf(",\n");
        }
    }
    flint_printf("]\n");
}
