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
acb_theta_agm_ctx_init_ext(acb_theta_agm_ctx_t ctx, acb_srcptr z,
                           const acb_mat_t tau)
{
    slong g = acb_mat_nrows(tau);
    slong nb = 1 << g;
    slong dim = 2 * (nb - 1);

    acb_theta_agm_ctx_is_ext(ctx) = 1;
    acb_theta_agm_ctx_dim(ctx) = dim;
    acb_theta_agm_ctx_init_internal(ctx, nb, g);

    acb_mat_set(acb_theta_agm_ctx_tau(ctx), tau);
    _acb_vec_set(acb_theta_agm_ctx_z(ctx), z, g);
}
