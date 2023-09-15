/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void acb_theta_g2_slash_basic_covariants(acb_poly_struct* res, const acb_mat_t c,
    const acb_poly_struct* cov, slong prec)
{
    slong klist[] = ACB_THETA_G2_BASIC_K;
    slong jlist[] = ACB_THETA_G2_BASIC_J;
    slong nb = ACB_THETA_G2_BASIC_NB;
    slong i;

    for (i = 0; i < nb; i++)
    {
        acb_theta_g2_detk_symj(&res[i], c, &cov[i], klist[i], jlist[i], prec);
    }
}


