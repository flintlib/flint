/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "mpn_extras.h"
#include "flint-impl.h"

void flint_mpn_debug(ulong_srcptr x, mp_size_t xsize)
{
    int i, j;
    char byte[9];
    byte[8] = 0;

    printf("\n");
    for (i = 0; i < xsize; i++)
    {
        printf("DIGIT %3d/" WORD_FMT "d: ", i, xsize);
        for (j = 0; j < FLINT_BITS; j++)
        {
            byte[j % 8] = (x[i] & (UWORD(1)<<j)) ? '1' : '0';
            if (j % 8 == 7)
                printf("%s ", byte);
        }
        printf(" (" WORD_FMT "u)\n", x[i]);
    }
}
