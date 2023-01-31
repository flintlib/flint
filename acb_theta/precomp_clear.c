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
acb_theta_precomp_clear(acb_theta_precomp_t D)
{
    slong g = acb_mat_nrows(acb_theta_precomp_exp_mat(D));
    slong nb_pow = D->indices[g];

    acb_mat_clear(acb_theta_precomp_exp_mat(D));
    flint_free(D->indices);
    if (nb_pow > 0)
        _acb_vec_clear(D->sqr_powers, nb_pow);
    _acb_vec_clear(D->exp_z, g * acb_theta_precomp_nb_z(D));
}
