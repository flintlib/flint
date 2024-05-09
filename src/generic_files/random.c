/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"

void _flint_rand_init_gmp_state(flint_rand_t state)
{
    FLINT_ASSERT(state->__gmp_state == NULL);
    state->__gmp_state = flint_malloc(sizeof(__gmp_randstate_struct));
    gmp_randinit_default(state->__gmp_state);
}

void _flint_rand_clear_gmp_state(flint_rand_t state)
{
    FLINT_ASSERT(state->__gmp_state != NULL);
    gmp_randclear(state->__gmp_state);
    flint_free(state->__gmp_state);
}

/* NOTE: The following functions are deprecated, and only exists for
   backwards-compatibility. */
void flint_randinit(flint_rand_t state)
{
    flint_rand_init(state);
}

void flint_randclear(flint_rand_t state)
{
    flint_rand_clear(state);
}

void flint_randseed(flint_rand_t state, ulong s0, ulong s1)
{
    flint_rand_set_seed(state, s0, s1);
}

void flint_get_randseed(ulong * s0, ulong * s1, flint_rand_t state)
{
    flint_rand_get_seed(s0, s1, state);
}

void _flint_rand_init_gmp(flint_rand_t state)
{
    if (state->__gmp_state == NULL)
        _flint_rand_init_gmp_state(state);
}

flint_rand_struct * flint_rand_alloc(void)
{
    return (flint_rand_struct *) flint_malloc(sizeof(flint_rand_struct));
}

void flint_rand_free(flint_rand_struct * state)
{
    flint_free(state);
}
