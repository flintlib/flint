/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_elliptic.h"

void
acb_elliptic_k(acb_t k, const acb_t m, slong prec)
{
    acb_t t;
    acb_init(t);
    acb_sub_ui(t, m, 1, prec);
    acb_neg(t, t);
    acb_sqrt(t, t, prec);
    acb_agm1(k, t, prec);
    acb_const_pi(t, prec);
    acb_div(k, t, k, prec);
    acb_mul_2exp_si(k, k, -1);
    acb_clear(t);
}

