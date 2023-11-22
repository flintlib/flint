/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"
#include "acb_theta.h"

#define ACB_THETA_G2_COV_J {6,0,4,8,2,6,8,12,0,4,6,10,2,4,8,0,6,6,2,4,2,4,0,2,2,0}

TEST_FUNCTION_START(acb_theta_g2_covariants_lead, state)
{
    slong iter;

    /* Test: agrees with g2_covariants */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong prec = 200 + n_randint(state, 100);
        slong bits = 2;
        slong nb = ACB_THETA_G2_COV_NB;
        slong jlist[] = ACB_THETA_G2_COV_J;
        acb_poly_struct * cov;
        acb_ptr res, test;
        acb_poly_t r;
        slong k;

        cov = flint_malloc(nb * sizeof(acb_poly_struct));
        for (k = 0; k < nb; k++)
        {
            acb_poly_init(&cov[k]);
        }
        acb_poly_init(r);
        res = _acb_vec_init(nb);
        test = _acb_vec_init(nb);

        acb_poly_randtest(r, state, 7, prec, bits);

        acb_theta_g2_covariants(cov, r, prec);
        acb_theta_g2_covariants_lead(res, r, prec);
        for (k = 0; k < nb; k++)
        {
            acb_poly_get_coeff_acb(&test[k], &cov[k], jlist[k]);
        }

        if (!_acb_vec_overlaps(res, test, nb))
        {
            flint_printf("FAIL\n");
            _acb_vec_printd(res, nb, 5);
            _acb_vec_printd(test, nb, 5);
            flint_printf("\ncovariants:\n");
            for (k = 0; k < nb; k++)
            {
                acb_poly_printd(&cov[k], 5);
                flint_printf("\n");
            }
            flint_abort();
        }

        for (k = 0; k < nb; k++)
        {
            acb_poly_clear(&cov[k]);
        }
        flint_free(cov);
        acb_poly_clear(r);
        _acb_vec_clear(res, nb);
        _acb_vec_clear(test, nb);
    }

    TEST_FUNCTION_END(state);
}
