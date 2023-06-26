/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

slong
acb_theta_balance_lowprec(slong g, slong prec)
{
    arb_t x;
    slong lowprec;

    arb_init(x);

    arb_set_si(x, prec);
    arb_root_ui(x, x, g + 1, prec);
    arb_mul_si(x, x, ACB_THETA_BALANCE_LOWPREC_MUL, prec);
    lowprec = arf_get_si(arb_midref(x), ARF_RND_CEIL);

    arb_clear(x);
    return lowprec;
}
