/*
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_factor.h"

/* a = (b + c) mod n

    Not a normal add mod function, assumes n is normalized (highest bit set)
    and b and c are reduced modulo n

*/

void
fmpz_factor_ecm_addmod(mp_ptr a, mp_ptr b, mp_ptr c, mp_ptr n, mp_limb_t n_size)
{
	mp_limb_t cy;

    cy = mpn_add_n(a, b, c, n_size);

    if (cy || (mpn_cmp(a, n, n_size) > 0))
        mpn_sub_n(a, a, n, n_size);
}
