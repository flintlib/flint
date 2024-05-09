/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpfr_vec.h"

void
_mpfr_vec_randtest(mpfr_ptr f, flint_rand_t state, slong len)
{
    slong i;

    if (!FLINT_RAND_GMP_STATE_IS_INITIALISED(state))
        _flint_rand_init_gmp_state(state);

    for (i = 0; i < len; i++)
        mpfr_urandomb(f + i, state->__gmp_state);
}
