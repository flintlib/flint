/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "ulong_extras.h"
#include "acb_mat.h"
#include "acb_theta.h"

int main(void)
{
    slong prec;
    acb_mat_t tau;
    acb_ptr dth, z;
    slong nb = 3 * 16;

    acb_mat_init(tau, 2, 2);
    dth = _acb_vec_init(nb);
    z = _acb_vec_init(2);

    acb_onei(acb_mat_entry(tau, 0, 0));
    acb_onei(acb_mat_entry(tau, 1, 1));
    acb_onei(acb_mat_entry(tau, 1, 0));
    acb_mul_2exp_si(acb_mat_entry(tau, 1, 0), acb_mat_entry(tau, 1, 0), -2);
    acb_set(acb_mat_entry(tau, 0, 1), acb_mat_entry(tau, 1, 0));
    acb_mat_printd(tau, 5);

    for (prec = 32; prec <= n_pow(2, 15); prec *= 2)
    {
        flint_printf("prec = %wd, naive:\n", prec);

        TIMEIT_START;
        acb_theta_jet_naive_all(dth, z, tau, 1, prec);
        TIMEIT_STOP;

        acb_printd(&dth[0], 5);
        flint_printf("\n");
        flint_printf("prec = %wd, ql:\n", prec);

        TIMEIT_START;
        acb_theta_jet_all(dth, z, tau, 1, prec);
        TIMEIT_STOP;

        acb_printd(&dth[0], 5);
        flint_printf("\n\n");
    }

    acb_mat_clear(tau);
    _acb_vec_clear(dth, nb);
    _acb_vec_clear(z, 2);
}
