/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "clock.h"
#include "tune.h"

ulong n_xgcd_0(nn_ptr, nn_ptr, ulong, ulong);
ulong n_xgcd_1(nn_ptr, nn_ptr, ulong, ulong);

#define DEFINE_IT(name, func)                   \
double name(void * vparam)                      \
{                                               \
    struct n_param_0 * param = vparam;          \
    nn_ptr ap, bp, xp, yp;                      \
    slong len;                                  \
    flint_time_t t0, t1;                        \
    slong ix;                                   \
                                                \
    ap = param->ap;                             \
    bp = param->bp;                             \
    xp = param->xp;                             \
    yp = param->yp;                             \
    len = param->len;                           \
                                                \
    flint_time_get(t0);                         \
    for (ix = 0; ix < len; ix++)                \
        func(ap, bp, xp[ix], yp[ix]);           \
    flint_time_get(t1);                         \
                                                \
    return flint_time_nsec_diff(t1, t0);        \
}

DEFINE_IT(_tune_n_xgcd_0, n_xgcd_0)
DEFINE_IT(_tune_n_xgcd_1, n_xgcd_0)
