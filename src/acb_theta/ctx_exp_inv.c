/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_theta.h"

void acb_theta_ctx_exp_inv(acb_t exp_inv, const acb_t exp, const acb_t z,
    int z_is_real, slong prec)
{
    if (z_is_real)
    {
        acb_conj(exp_inv, exp);
    }
    else if (acb_contains_zero(exp)
        || acb_rel_error_bits(exp) <= 0.75 * prec)
    {
        acb_t x;
        acb_init(x);
        acb_neg(x, z);
        acb_exp_pi_i(exp_inv, x, prec);
        acb_clear(x);
    }
    else
    {
        acb_inv(exp_inv, exp, prec);
    }
}
