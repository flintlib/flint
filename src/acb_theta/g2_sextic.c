/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_theta.h"

void
acb_theta_g2_sextic(acb_poly_t res, const acb_mat_t tau, slong prec)
{
    acb_t chi5;

    acb_init(chi5);
    acb_theta_g2_sextic_chi5(res, chi5, tau, prec);
    acb_clear(chi5);
}
