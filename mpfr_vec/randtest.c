/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpfr_vec.h"

void
_mpfr_vec_randtest(flint_mpfr * f, flint_rand_t state, slong len)
{
    slong i;

    _flint_rand_init_gmp(state);

    for (i = 0; i < len; i++)
        mpfr_urandomb(f + i, state->gmp_state);
}
