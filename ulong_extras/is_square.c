/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#define ulong ulongxx /* interferes with system includes */
#include <math.h>
#undef ulong
#define ulong mp_limb_t
#include <gmp.h>
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
   
int n_is_square(mp_limb_t x)
{
    mp_limb_t sq;

    if (!mod64[x % UWORD(64)]) return 0;
    if (!mod63[x % UWORD(63)]) return 0;
    if (!mod65[x % UWORD(65)]) return 0;

    sq = (mp_limb_t) (sqrt((double) x) + 0.5);
    
    return (x == sq*sq);
}
