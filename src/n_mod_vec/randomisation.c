/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_mod.h"
#include "n_mod_vec.h"

void _n_mod_vec_rand(nn_ptr rp, flint_rand_t state, slong len, n_mod_ctx_srcptr ctx)
{
    slong ix;

    for (ix = 0; ix < len; ix++)
        rp[ix] = n_mod_rand(state, ctx);
}
