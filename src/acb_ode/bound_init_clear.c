/*
    Copyright (C) 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_ode.h"
#include "acb_poly.h"
#include "arb_poly.h"

void
acb_ode_bound_init(acb_ode_bound_t bound)
{
    bound->prec = -1;
    arb_init(bound->cst);
    arb_poly_init(bound->den_lbound);
    bound->den_rt = NULL;
    bound->den_rt_len = 0;
    bound->all_nums = NULL;
    bound->all_nums_len = 0;
    bound->pol_part_len = 0;
}

void
acb_ode_bound_clear(acb_ode_bound_t bound)
{
    arb_clear(bound->cst);
    arb_poly_clear(bound->den_lbound);
    _arb_vec_clear(bound->den_rt, bound->den_rt_len);
    _acb_poly_vec_clear(bound->all_nums, bound->all_nums_len);
}
