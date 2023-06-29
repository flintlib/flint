/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

slong
acb_theta_naive_newprec(slong prec, slong coord, slong dist,
    slong max_dist, slong ord)
{
    double r = ((double) dist) / (max_dist + 2);
    double neg = r * r * prec;
    double pos = ord * n_clog(1 + FLINT_ABS(coord), 2);

    return FLINT_MAX(ACB_THETA_ELD_DEFAULT_PREC, ceil((double) prec - neg + pos));
}
