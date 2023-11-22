/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_eld_print(const acb_theta_eld_t E)
{
    slong d = acb_theta_eld_dim(E);
    slong g = acb_theta_eld_ambient_dim(E);
    slong k;

    for (k = 0; k < g - d; k++)
    {
        flint_printf("  ");
    }
    flint_printf("Slice (...");
    for (k = 0; k < g - d; k++)
    {
        flint_printf(", %wd", acb_theta_eld_coord(E, k + d));
    }
    flint_printf("): from %wd to %wd (mid: %wd)\n",
        acb_theta_eld_min(E), acb_theta_eld_max(E), acb_theta_eld_mid(E));
    if (d > 1)
    {
        for (k = 0; k < acb_theta_eld_nr(E); k++)
        {
            acb_theta_eld_print(acb_theta_eld_rchild(E, k));
        }
        for (k = 0; k < acb_theta_eld_nl(E); k++)
        {
            acb_theta_eld_print(acb_theta_eld_lchild(E, k));
        }
    }
}
