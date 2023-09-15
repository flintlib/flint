/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void acb_theta_g2_covariant_weight(slong* k, slong* j, const fmpz_mpoly_t pol,
    const fmpz_mpoly_ctx_t ctx)
{
    slong e[ACB_THETA_G2_BASIC_NB];
    slong klist[] = ACB_THETA_G2_BASIC_K;
    slong jlist[] = ACB_THETA_G2_BASIC_J;
    slong i;

    fmpz_mpoly_get_term_exp_si(e, pol, 0, ctx);
    *k = 0;
    *j = 0;
    for (i = 0; i < ACB_THETA_G2_BASIC_NB; i++)
    {
        *k += e[i] * klist[i];
        *j += e[i] * jlist[i];
    }
    *k -= (*j/2);
}
