/*
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "mpn_extras.h"

/* x = (a - b) mod n 

    Not a normal sub mod function, assumes n is normalized (higest bit set)
    and a and b are reduced modulo n

*/

void
fmpz_factor_ecm_submod(mp_ptr x, mp_ptr a, mp_ptr b, mp_ptr n, mp_limb_t n_size)
{
    mp_ptr temp;

    TMP_INIT;

    TMP_START;
    temp = TMP_ALLOC(n_size * sizeof(mp_limb_t));

    if (mpn_cmp(a, b, n_size) > 0)
        mpn_sub_n(x, a, b, n_size);
    else
    {
        mpn_sub_n(temp, n, b, n_size);
        mpn_add_n(x, temp, a, n_size);
    }

    TMP_END;
}
