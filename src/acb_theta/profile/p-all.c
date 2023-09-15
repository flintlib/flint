/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"
#include "profiler.h"

int main(void)
{
    slong prec;
    acb_mat_t tau;
    acb_ptr th, z;
    slong nb = 16;

    acb_mat_init(tau, 2, 2);
    th = _acb_vec_init(nb);
    z = _acb_vec_init(2);

    acb_onei(acb_mat_entry(tau, 0, 0));
    acb_onei(acb_mat_entry(tau, 1, 1));
    acb_onei(acb_mat_entry(tau, 1, 0));
    acb_mul_2exp_si(acb_mat_entry(tau, 1, 0), acb_mat_entry(tau, 1, 0), -2);
    acb_set(acb_mat_entry(tau, 0, 1), acb_mat_entry(tau, 1, 0));
    acb_mat_printd(tau, 5);

    for (prec = 32; prec <= n_pow(2, 14); prec *= 2)
    {
        flint_printf("prec = %wd, naive:\n", prec);

        TIMEIT_START

            acb_theta_naive_all(th, z, 1, tau, prec);

        TIMEIT_STOP;

        acb_printd(&th[0], 5);
        flint_printf("\n");
        flint_printf("prec = %wd, ql:\n", prec);

        TIMEIT_START

            acb_theta_all(th, z, tau, 0, prec);

        TIMEIT_STOP;
        acb_printd(&th[0], 5);
        flint_printf("\n\n");
    }

    acb_mat_clear(tau);
    _acb_vec_clear(th, nb);
    _acb_vec_clear(z, 2);
}
