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
acb_theta_char_dot_acb(acb_t x, ulong a, acb_srcptr z, slong g, slong prec)
{
    slong* v;
    slong k;

    v = flint_malloc(g * sizeof(slong));
    
    acb_theta_char_get_slong(v, a, g);
    acb_dot_si(x, NULL, 0, z, 1, v, 1, g, prec);
    
    flint_free(v);
}
