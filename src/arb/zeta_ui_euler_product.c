/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "acb_dirichlet.h"

void
arb_zeta_inv_ui_euler_product(arb_t z, ulong s, slong prec)
{
    const signed char chi[1] = {1};
    _acb_dirichlet_euler_product_real_ui(z, s, chi, 1, 1, prec);
}

void
arb_zeta_ui_euler_product(arb_t z, ulong s, slong prec)
{
    const signed char chi[1] = {1};
    _acb_dirichlet_euler_product_real_ui(z, s, chi, 1, 0, prec);
}

void
arb_zeta_ui_asymp(arb_t z, ulong s, slong prec)
{
    const signed char chi[1] = {1};
    _acb_dirichlet_euler_product_real_ui(z, s, chi, 1, 0, prec);
}

