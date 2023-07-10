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
acb_theta_get_a0(acb_ptr r, acb_srcptr th, slong g)
{    
    acb_ptr v;
    slong a;
    slong n = 1 << g;

    v = _acb_vec_init(n);    
    for (a = 0; a < n; a++)
    {
        acb_set(&v[a], &th[n * a]);
    }
    _acb_vec_set(r, v, n);
    _acb_vec_clear(v, n);
}
