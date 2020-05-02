/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

void
ca_qqbar_cot_pi(ca_qqbar_t res, slong p, ulong q)
{
    ca_qqbar_t t, u;

    ca_qqbar_init(t);
    ca_qqbar_init(u);

    ca_qqbar_sin_pi(t, p, q);
    ca_qqbar_cos_pi(u, p, q);

    /* the monic polynomials have smaller coefficients */
    ca_qqbar_mul_2exp_si(t, t, 1);
    ca_qqbar_mul_2exp_si(u, u, 1);

    ca_qqbar_div(res, u, t);

    ca_qqbar_clear(t);
    ca_qqbar_clear(u);
}

