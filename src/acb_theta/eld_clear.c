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
acb_theta_eld_clear(acb_theta_eld_t E)
{
    slong k;
    slong nr = acb_theta_eld_nr(E);
    slong nl = acb_theta_eld_nl(E);

    if (nr > 0)
    {
        for (k = 0; k < nr; k++)
        {
            acb_theta_eld_clear(acb_theta_eld_rchild(E, k));
        }
        flint_free(E->rchildren);
    }
    if (nl > 0)
    {
        for (k = 0; k < nl; k++)
        {
            acb_theta_eld_clear(acb_theta_eld_lchild(E, k));
        }
        flint_free(E->lchildren);
    }

    flint_free(E->last_coords);
    flint_free(E->box);
}
