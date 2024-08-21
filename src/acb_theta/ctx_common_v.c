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

static void
_arb_vec_union(arb_ptr res, arb_srcptr v1, arb_srcptr v2, slong len, slong prec)
{
    slong j;

    for (j = 0; j < len; j++)
    {
        arb_union(&res[j], &v1[j], &v2[j], prec);
    }
}

void acb_theta_ctx_common_v(arb_ptr v, const acb_theta_ctx_t ctx, slong prec)
{
    slong g = acb_theta_ctx_g(ctx);
    slong nb = acb_theta_ctx_nb(ctx);
    slong j;

    if (nb == 0) return;

    _arb_vec_set(v, acb_theta_ctx_vs(ctx), g);
    for (j = 1; j < nb; j++)
    {
	_arb_vec_union(v, v, acb_theta_ctx_vs(ctx) + j * g, g, prec);
    }
}

