/*
    Copyright (C) 2025 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "arb.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_reduce_z, state)
{
    slong iter;

    /* Test: entries of r are even integers, im(new_zs) is im(zs - tau * r) */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong nb = n_randint(state, 3);
        slong prec = 100 + n_randint(state, 200);
        slong bits = n_randint(state, 4);
        acb_ptr zs, new_zs, cs;
        acb_mat_t tau;
        arb_mat_t im;
        arb_ptr rs, test, y;
        fmpz_t x;
        slong j;
        int res;

        acb_mat_init(tau, g, g);
        arb_mat_init(im, g, g);
        zs = _acb_vec_init(nb * g);
        new_zs = _acb_vec_init(nb * g);
        cs = _acb_vec_init(nb);
        test = _arb_vec_init(g);
        y = _arb_vec_init(g);
        rs = _arb_vec_init(nb * g);
        fmpz_init(x);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_siegel_randtest_vec(zs, state, nb * g, prec);
        acb_mat_get_imag(im, tau);

        res = acb_theta_reduce_z(new_zs, rs, cs, zs, nb, tau, prec);

        if (res)
        {
            for (j = 0; j < g * nb; j++)
            {
                res = arb_get_unique_fmpz(x, &rs[j]);
                if (!res || !arb_is_exact(&rs[j]) || fmpz_mod_ui(x, x, 2) != 0)
                {
                    flint_printf("FAIL (integers)\n");
                    flint_abort();
                }
            }

            for (j = 0; j < nb; j++)
            {
                arb_mat_vector_mul_col(test, im, rs + j * g, prec);
                _acb_vec_get_imag(y, zs + j * g, g);
                _arb_vec_sub(test, y, test, g, prec);
                _acb_vec_get_imag(y, new_zs + j * g, g);
                if (!_arb_vec_overlaps(test, y, g))
                {
                    flint_printf("FAIL (overlap)\n");
                    _arb_vec_printd(y, g, 5);
                    _arb_vec_printd(test, g, 5);
                    flint_abort();
                }
            }
        }

        acb_mat_clear(tau);
        arb_mat_clear(im);
        _acb_vec_clear(zs, nb * g);
        _acb_vec_clear(new_zs, nb * g);
        _acb_vec_clear(cs, nb);
        _arb_vec_clear(test, g);
        _arb_vec_clear(y, g);
        _arb_vec_clear(rs, nb * g);
        fmpz_clear(x);
    }

    TEST_FUNCTION_END(state);
}
