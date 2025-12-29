/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_poly.h"
#include "acb_dirichlet.h"

void
acb_hurwitz_zeta(acb_t z, const acb_t s, const acb_t a, slong prec)
{
    acb_dirichlet_hurwitz(z, s, a, prec);
}

void
acb_zeta(acb_t z, const acb_t s, slong prec)
{
    acb_dirichlet_zeta(z, s, prec);
}
