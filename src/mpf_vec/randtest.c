/*
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpf_vec.h"

void
_mpf_vec_randtest(mpf * f, flint_rand_t state, slong len, flint_bitcnt_t bits)
{
    slong i;

    _flint_rand_init_gmp(state);

    for (i = 0; i < len; i++)
        mpf_urandomb(f + i, state->gmp_state, bits);
}
