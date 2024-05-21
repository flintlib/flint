/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "clock.h"
#include "tune.h"
#include "n_mod.h"
#include "n_mod_vec.h"

void _n_mod_vec_add_0(nn_ptr, nn_srcptr, nn_srcptr, slong, ulong);
void _n_mod_vec_sub_0(nn_ptr, nn_srcptr, nn_srcptr, slong, ulong);
void _n_mod_vec_add_1(nn_ptr, nn_srcptr, nn_srcptr, slong, ulong);
void _n_mod_vec_sub_1(nn_ptr, nn_srcptr, nn_srcptr, slong, ulong);

#define DEFINE_IT(name, func)                   \
double name(void * vparam)                      \
{                                               \
    struct n_mod_vec_param_0 * param = vparam;  \
    nn_ptr rp, ap, bp;                          \
    slong len;                                  \
    ulong nu;                                   \
    flint_time_t t0, t1;                        \
                                                \
    rp = param->rp;                             \
    ap = param->ap;                             \
    bp = param->bp;                             \
    len = param->len;                           \
    nu = param->ctx->nu;                        \
                                                \
    flint_time_get(t0);                         \
    func(rp, ap, bp, len, nu);                  \
    flint_time_get(t1);                         \
                                                \
    return flint_time_nsec_diff(t1, t0);        \
}

DEFINE_IT(_tune_n_mod_vec_add_0, _n_mod_vec_add_0)
DEFINE_IT(_tune_n_mod_vec_add_1, _n_mod_vec_add_1)
DEFINE_IT(_tune_n_mod_vec_sub_0, _n_mod_vec_sub_0)
DEFINE_IT(_tune_n_mod_vec_sub_1, _n_mod_vec_sub_1)
