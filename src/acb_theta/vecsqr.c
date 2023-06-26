/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_vecsqr(acb_ptr th2, acb_srcptr th, slong n, slong prec)
{
    slong k;
    for (k = 0; k < n; k++)
    {
        acb_sqr(&th2[k], &th[k], prec);
    }
}
