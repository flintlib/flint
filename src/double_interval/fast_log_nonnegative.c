/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "double_interval.h"
#include "mag.h"

double mag_d_log_lower_bound(double x);
double mag_d_log_upper_bound(double x);

di_t di_fast_log_nonnegative(di_t x)
{
    di_t res;

    if (x.a <= 0.0)
        res.a = -D_INF;
    else
        res.a = mag_d_log_lower_bound(x.a);

    res.b = mag_d_log_upper_bound(x.b);

    return res;
}
