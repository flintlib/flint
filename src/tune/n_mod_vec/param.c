/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "tune.h"
#include "n_mod.h"
#include "n_mod_vec.h"

#if FLINT64
# define N_0 UWORD(7365182178263871635)
#else
# define N_0 UWORD(1236571635)
#endif

void * n_mod_vec_param_init_generate_0(void)
{
    struct n_mod_vec_param_0 * param;
    nn_ptr rp, ap, bp;
    slong len;
    ulong n;

    param = flint_malloc(sizeof(struct n_mod_vec_param_0));

    len = 1000;
    rp = flint_malloc(sizeof(ulong) * len);
    ap = flint_malloc(sizeof(ulong) * len);
    bp = flint_malloc(sizeof(ulong) * len);
    n = N_0;

    param->rp = rp;
    param->ap = ap;
    param->bp = bp;
    param->len = len;
    n_mod_ctx_init(param->ctx, n);

    return param;
}

void n_mod_vec_param_clear(void * vparam)
{
    struct n_mod_vec_param_0 * param = vparam;

    flint_free(param->rp);
    flint_free(param->ap);
    flint_free(param->bp);
    n_mod_ctx_clear(param->ctx);
    flint_free(param);
}
