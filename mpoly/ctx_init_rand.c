/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

void mpoly_ctx_init_rand(mpoly_ctx_t mctx, flint_rand_t state, slong max_nvars)
{
    ordering_t ord;
    slong nvars;

    ord = mpoly_ordering_randtest(state);
    nvars = n_randint(state, max_nvars + 1);
    mpoly_ctx_init(mctx, nvars, ord);
}
