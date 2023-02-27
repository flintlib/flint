/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2009, 2010 Andy Novocin
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_lll.h"

int
fmpz_lll_shift(const fmpz_mat_t B)
{
    int i, n = B->c;
    int shift = 0;
    for (i = 0; i < B->r; i++)
    {
        int j;
        for (j = n - 1;
             j >= i + shift + 1 && fmpz_size(fmpz_mat_entry(B, i, j)) == 0L;
             j--) ;

        if (shift < j - i)
            shift = j - i;

    }

    return shift;
}
