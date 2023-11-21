/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_hurwitz_precomp_clear(acb_dirichlet_hurwitz_precomp_t pre)
{
    acb_clear(&pre->s);

    if (pre->A != 0)
    {
        mag_clear(&pre->err);
        _acb_vec_clear(pre->coeffs, pre->N * pre->K);
    }
}

