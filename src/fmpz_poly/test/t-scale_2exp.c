/*
    Copyright (C) 2016, Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"
TEST_FUNCTION_START(fmpz_poly_scale_2exp, state)
{
    int iter;

    /* Check that scale(k) followed by scale(-k) is identity */
    for (iter = 0; iter < 1000; iter++)
    {
        fmpz_poly_t f, g;
        slong k;

        k = n_randint(state, 200);

        fmpz_poly_init(f);
        fmpz_poly_init(g);

        fmpz_poly_randtest(f, state, n_randint(state, 100), 200);
        _fmpz_poly_remove_content_2exp(f->coeffs, f->length);
        fmpz_poly_set(g, f);

        _fmpz_poly_scale_2exp(f->coeffs, f->length, k);
        _fmpz_poly_scale_2exp(f->coeffs, f->length, -k);

        if ( !fmpz_poly_equal(f, g) )
        {
            flint_printf("FAIL:\n");
            fmpz_poly_print(f); flint_printf("\n\n");
            flint_printf("ERROR\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
    }

    /* Check consistency with polynomial evaluation */
    for (iter = 0; iter < 1000; iter++)
    {
        slong k;
        fmpz_poly_t f, g;
        fmpz_t a1,a2,res1,res2;

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_init(a1);
        fmpz_init(a2);
        fmpz_init(res1);
        fmpz_init(res2);

        k = n_randint(state, 30);

        fmpz_randtest(a2, state, 100);
        fmpz_mul_2exp(a1, a2, k);

        /* positive case k */
        fmpz_poly_set(g, f);
        _fmpz_poly_scale_2exp(g->coeffs, g->length, k);

        fmpz_poly_evaluate_fmpz(res1, f, a1);
        fmpz_fdiv_q_2exp(res1, res1, fmpz_val2(res1));

        fmpz_poly_evaluate_fmpz(res2, g, a2);
        fmpz_fdiv_q_2exp(res2, res2, fmpz_val2(res2));

        if ( !fmpz_equal(res1, res2) )
        {
            flint_printf("FAIL:\n");
            flint_printf("f = "); fmpz_poly_print(f); flint_printf("\n\n");
            flint_printf("g = "); fmpz_poly_print(g); flint_printf("\n\n");
            flint_printf("ERROR\n");
            fflush(stdout);
            flint_abort();
        }

        /* negative case -k */
        fmpz_poly_set(g, f);
        _fmpz_poly_scale_2exp(g->coeffs, g->length, -k);

        fmpz_poly_evaluate_fmpz(res1, f, a2);
        fmpz_fdiv_q_2exp(res2, res2, fmpz_val2(res2));

        fmpz_poly_evaluate_fmpz(res1, g, a1);
        fmpz_fdiv_q_2exp(res1, res1, fmpz_val2(res1));

        if ( !fmpz_equal(res1, res2) )
        {
            flint_printf("FAIL:\n");
            flint_printf("f = "); fmpz_poly_print(f); flint_printf("\n\n");
            flint_printf("g = "); fmpz_poly_print(g); flint_printf("\n\n");
            flint_printf("ERROR\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_clear(a1);
        fmpz_clear(a2);
        fmpz_clear(res1);
        fmpz_clear(res2);
    }

    TEST_FUNCTION_END(state);
}
