/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "flint.h"
#include "ulong_extras.h"

int mod64[64] = {1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,
                 0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,
                 0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0};

int mod65[65] = {1,1,0,0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,0,0,0,
                 0,0,0,1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1,0,0,0,
                 0,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,0,1,0,0,1};

int mod63[63] = {1,1,0,0,1,0,0,1,0,1,0,0,0,0,1,0,1,0,1,0,0,
                 0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,
                 0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0};

int n_is_square(ulong x)
{
    ulong sq;

    if (!mod64[x % UWORD(64)]) return 0;
    if (!mod63[x % UWORD(63)]) return 0;
    if (!mod65[x % UWORD(65)]) return 0;

    sq = (ulong) (sqrt((double) x) + 0.5);

    return (x == sq*sq);
}
