/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "tune.h"

void * n_param_init_generate_0(void)
{
    struct n_param_0 * param;
    nn_ptr ap, bp, xp, yp;
    slong len;
    slong ix;
    flint_rand_t state;

    param = flint_malloc(sizeof(struct n_param_0));
    flint_rand_init(state);

    len = 100;
    ap = flint_malloc(sizeof(ulong) * len);
    bp = flint_malloc(sizeof(ulong) * len);
    xp = flint_malloc(sizeof(ulong) * len);
    yp = flint_malloc(sizeof(ulong) * len);

    for (ix = 0; ix < len; ix++)
    {
        xp[ix] = n_randlimb(state);
        yp[ix] = n_randlimb(state);
        if (xp[ix] < yp[ix])
            FLINT_SWAP(ulong, xp[ix], yp[ix]);
    }

    param->ap = ap;
    param->bp = bp;
    param->xp = xp;
    param->yp = yp;
    param->len = len;

    flint_rand_clear(state);

    return param;
}

void n_param_clear(void * vparam)
{
    struct n_param_0 * param = vparam;

    flint_free(param->ap);
    flint_free(param->bp);
    flint_free(param->xp);
    flint_free(param->yp);
    flint_free(param);
}
