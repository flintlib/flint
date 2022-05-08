/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "gmp.h"

void _flint_rand_init_gmp(flint_rand_t state)
{
    if (!state->gmp_init)
    {
        gmp_randinit_default((__gmp_randstate_struct *) state->gmp_state);
        state->gmp_init = 1;
    }
}

void flint_randclear(flint_rand_t state)
{
    if (state->gmp_init)
        gmp_randclear((__gmp_randstate_struct *) state->gmp_state);
}
