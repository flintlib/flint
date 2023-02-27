/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2009 Thomas Boothby

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

#define FLINT_ONE_LINE_MULTIPLIER 480

mp_limb_t n_factor_one_line(mp_limb_t n, ulong iters)
{
    mp_limb_t orig_n = n, in, square, sqrti, mod, factor, factoring = iters, iin;
    n *= FLINT_ONE_LINE_MULTIPLIER;

    iin = 0;
    in = n;
    while (factoring && (iin < in))
    {
        sqrti = n_sqrt(in);
        sqrti++;
        square = sqrti*sqrti;
        mod = square - in;
        if (n_is_square(mod)) 
        {
            factor = n_sqrt(mod);
            sqrti -= factor;
            factor = n_gcd(orig_n, sqrti);
            if (factor != UWORD(1)) 
            { 
                return factor;
            }
        }     
        factoring--;    
        iin = in;
        in += n;
    }

    return 0;
}
