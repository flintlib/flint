/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_agm_ctx_reset_steps(acb_theta_agm_ctx_t ctx, slong k, slong m)
{
    slong g = acb_theta_agm_ctx_g(ctx);
    slong prev = acb_theta_agm_ctx_nb_bad(ctx, k);
    slong nb_th;

    nb_th = 1 << g;
    if (acb_theta_agm_ctx_is_ext(ctx))
        nb_th *= 2;

    if (prev == 0 && m > 0)
    {
        acb_theta_agm_ctx_nb_bad(ctx, k) = m;
        acb_theta_agm_ctx_roots(ctx, k) = _acb_vec_init(m * nb_th);
    }
    else if (prev > 0 && m == 0)
    {
        acb_theta_agm_ctx_nb_bad(ctx, k) = 0;
        _acb_vec_clear(acb_theta_agm_ctx_roots(ctx, k), prev * nb_th);
    }
    else if (prev > 0 && m > 0)
    {
        acb_theta_agm_ctx_reset_steps(ctx, k, 0);
        acb_theta_agm_ctx_reset_steps(ctx, k, m);
    }
}
