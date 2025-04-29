/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_eld_clear(acb_theta_eld_t E)
{
    slong k;

    if ((E->nr) > 0)
    {
        for (k = 0; k < (E->nr); k++)
        {
            acb_theta_eld_clear(&E->rchildren[k]);
        }
        flint_free(E->rchildren);
    }
    if ((E->nl) > 0)
    {
        for (k = 0; k < (E->nl); k++)
        {
            acb_theta_eld_clear(&E->lchildren[k]);
        }
        flint_free(E->lchildren);
    }

    flint_free(E->last_coords);
    flint_free(E->box);
}
