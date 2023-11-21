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
#include "acb_modular.h"

void
acb_elliptic_invariants(acb_t g2, acb_t g3, const acb_t tau, slong prec)
{
    acb_struct t[2];

    acb_init(t);
    acb_init(t + 1);

    acb_modular_eisenstein(t, tau, 2, prec);

    acb_mul_ui(g2, t, 60, prec);
    acb_mul_ui(g3, t + 1, 140, prec);

    acb_clear(t);
    acb_clear(t + 1);
}

