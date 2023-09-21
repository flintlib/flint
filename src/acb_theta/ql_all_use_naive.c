/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int acb_theta_ql_all_use_naive(slong g, slong prec)
{
    if (g <= 2)
    {
        return (prec <= 1500);
    }
    else
    {
        return 0;
    }
}