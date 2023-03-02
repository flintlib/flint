/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "mpn_extras.h"


void flint_mpn_debug(mp_srcptr x, mp_size_t xsize)
{
    int i, j;
    char byte[9];
    byte[8] = 0;

    flint_printf("\n");
    for (i = 0; i < xsize; i++)
    {
        flint_printf("DIGIT %3d/%wd: ", i, xsize);
        for (j = 0; j < FLINT_BITS; j++)
        {
            byte[j % 8] = (x[i] & (UWORD(1)<<j)) ? '1' : '0';
            if (j % 8 == 7)
                flint_printf("%s ", byte);
        }
        flint_printf(" (%wu)\n", x[i]);
    }
}
