/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"
#include "acb_modular.h"
#include "acb_theta.h"

/* TODO: work in progress */

void
acb_theta_sum_fixed_ab(acb_ptr th, ulong ab, const acb_theta_ctx_t ctx, slong prec)
{
    slong g = acb_theta_ctx_g(ctx);
    slong nb = acb_theta_ctx_nb(ctx);
    slong j;

    if (nb == 0)
    {
	return;
    }

    if (g == 1)
    {
	acb_ptr res;

	res = _acb_vec_init(4);
	for (j = 0; j < nb; j++)
	{
	    /* acb_modular_theta_sum recomputes the inverse of exp_z */
	    /* todo: store w_is_unit as part of context */
	    acb_modular_theta_sum(&res[0], &res[1], &res[2], &res[3],
		&acb_theta_ctx_exp_zs(ctx)[j], 0,
		acb_mat_entry(acb_theta_ctx_exp_tau(ctx), 0, 0), 1, prec);
	    if (ab <= 1)
	    {
		acb_set(&th[j], &res[2 + ab]);
	    }
	    else
	    {
		if (ab == 2)
		{
		    acb_set(&th[j], &res[1]);
		}
		else
		{
		    acb_neg(&th[j], &res[0]);
		}
		acb_mul(&th[j], &th[j], acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), 0, 0), prec);
	    }
	}
	_acb_vec_clear(res, 4);
    }
    else
    {
	
    }
}
