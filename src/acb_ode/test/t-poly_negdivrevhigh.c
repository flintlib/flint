/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Generated using Claude Opus 4.8 */

#include "test_helpers.h"
#include "acb.h"
#include "acb_ode.h"
#include "acb_types.h"

TEST_FUNCTION_START(acb_ode_poly_negdivrevhigh, state)
{
    for (slong iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong len = n_randint(state, 12);
        slong prec = 2 + n_randint(state, 200);
        slong bits = 1 + n_randint(state, 8);
        int use_cst = n_randint(state, 4);

        acb_ptr a = _acb_vec_init(len);
        acb_ptr b = _acb_vec_init(len);
        acb_ptr res = _acb_vec_init(len);
        acb_t cst, acc;

        acb_init(cst);
        acb_init(acc);

        for (slong i = 0; i < len; i++)
        {
            acb_randtest(a + i, state, prec, bits);
            acb_randtest(b + i, state, prec, bits);
        }

        if (use_cst)
            acb_randtest(cst, state, prec, bits);

        _acb_ode_poly_negdivrevhigh(res, a, use_cst ? cst : NULL, b, len, prec);

        /* The coefficient of X^d in a(1/X)*res(X) is sum_i a[i]*res[d+i], so
           we must have sum_i a[i]*res[d+i] + cst*b[d] = 0 for all d. */
        for (slong d = 0; d < len; d++)
        {
            if (use_cst)
                acb_mul(acc, cst, b + d, prec);
            else
                acb_set(acc, b + d);

            acb_dot(acc, acc, 0, a, 1, res + d, 1, len - d, prec);

            if (!acb_contains_zero(acc))
            {
                flint_printf("FAIL (identity, d = %wd)\n\n", d);
                flint_printf("len = %wd, use_cst = %d\n", len, use_cst);
                flint_printf("a = %{acb*}\n", a, len);
                flint_printf("b = %{acb*}\n", b, len);
                if (use_cst)
                    flint_printf("cst = %{acb}\n", cst);
                flint_printf("res = %{acb*}\n", res, len);
                flint_printf("residue = %{acb}\n", acc);
                flint_abort();
            }
        }

        acb_clear(acc);
        acb_clear(cst);
        _acb_vec_clear(res, len);
        _acb_vec_clear(b, len);
        _acb_vec_clear(a, len);
    }

    TEST_FUNCTION_END(state);
}
