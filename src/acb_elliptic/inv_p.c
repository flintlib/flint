/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_elliptic.h"

void
acb_elliptic_inv_p(acb_t res, const acb_t z, const acb_t tau, slong prec)
{
    acb_t e1, e2, e3;

    acb_init(e1);
    acb_init(e2);
    acb_init(e3);

    acb_elliptic_roots(e1, e2, e3, tau, prec);

    acb_sub(e1, z, e1, prec);
    acb_sub(e2, z, e2, prec);
    acb_sub(e3, z, e3, prec);

    acb_elliptic_rf(res, e1, e2, e3, 0, prec);

    acb_clear(e1);
    acb_clear(e2);
    acb_clear(e3);
}

