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
acb_dirichlet_gauss_sum_order2(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
{
    if (dirichlet_parity_char(G, chi))
    {
        arb_zero(acb_realref(res));
        arb_sqrt_ui(acb_imagref(res), G->q, prec);
    }
    else
    {
        arb_zero(acb_imagref(res));
        arb_sqrt_ui(acb_realref(res), G->q, prec);
    }
}
