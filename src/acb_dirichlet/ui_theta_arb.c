/*
    Copyright (C) 2016 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_ui_theta_arb(acb_t res, const dirichlet_group_t G, ulong a, const arb_t t, slong prec)
{
    dirichlet_char_t chi;

    dirichlet_char_init(chi, G);
    dirichlet_char_log(chi, G, a);

    acb_dirichlet_theta_arb(res, G, chi, t, prec);

    dirichlet_char_clear(chi);
}
