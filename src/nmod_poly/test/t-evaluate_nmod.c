/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_evaluate_nmod, state)
{
    int i, j, result = 1;

    /* Check evaluation at 1 gives sum of coeffs */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a;
        ulong n = n_randtest_not_zero(state);
        ulong sum, eval, eval_v1, eval_v2, eval_v3;

        nmod_poly_init(a, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));

        // main function
        eval = nmod_poly_evaluate_nmod(a, UWORD(1));

        // variants  (mulmod_precomp_shoup requires 1 < n)
        eval_v1 = _nmod_poly_evaluate_nmod(a->coeffs, a->length, UWORD(1), a->mod);
        if ((n > UWORD(1)) && (a->mod.norm > 0))
        {
            ulong one_precomp = n_mulmod_precomp_shoup(UWORD(1), n);
            eval_v2 = _nmod_poly_evaluate_nmod_precomp(a->coeffs, a->length, UWORD(1), one_precomp, a->mod);
#if FLINT_BITS == 64
            if (n <= UWORD(6148914691236517205))
#else // FLINT_BITS == 32
            if (n <= UWORD(1431655765))
#endif
            {
                eval_v3 = _nmod_poly_evaluate_nmod_precomp_lazy(a->coeffs, a->length, UWORD(1), one_precomp, a->mod);
            }
            else
            {
                eval_v3 = eval_v1;
            }
        }
        else
        {
            eval_v2 = eval_v1;
            eval_v3 = eval_v1;
        }

        sum = 0;
        for (j = 0; j < a->length; j++)
           sum = n_addmod(sum, nmod_poly_get_coeff_ui(a, j), n);

        result = (sum == eval) &&
                 (eval == eval_v1) && (eval == eval_v2) && (eval == eval_v3);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a->length = %wd, n = %wu\n", a->length, a->mod.n);
            flint_printf("sum = %wu, eval = %wu, eval_v1 = %wu, eval_v2 = %wu, eval_v3 = %wu\n",
                         sum, eval, eval_v1, eval_v2, eval_v3);
            nmod_poly_print(a), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
    }

    /* Check a(c) + b(c) = (a + b)(c) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, acopy;
        ulong n = n_randtest_not_zero(state);
        ulong eval1, eval2, c, eval1_v1, eval1_v2, eval1_v3, eval2_v1, eval2_v2, eval2_v3;

        nmod_poly_init(a, n);
        nmod_poly_init(acopy, n);
        nmod_poly_init(b, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        nmod_poly_randtest(b, state, n_randint(state, 100));
        nmod_poly_set(acopy, a);

        c = n_randint(state, n);

        // main function
        eval1 = nmod_poly_evaluate_nmod(a, c);
        eval1 = n_addmod(eval1, nmod_poly_evaluate_nmod(b, c), n);

        nmod_poly_add(a, a, b);
        eval2 = nmod_poly_evaluate_nmod(a, c);

        // variants
        nmod_poly_set(a, acopy);
        eval1_v1 = _nmod_poly_evaluate_nmod(a->coeffs, a->length, c, a->mod);
        eval1_v1 = n_addmod(eval1_v1, _nmod_poly_evaluate_nmod(b->coeffs, b->length, c, b->mod), n);

        if (a->mod.norm > 0)
        {
            ulong c_precomp = n_mulmod_precomp_shoup(c, n);
            eval1_v2 = _nmod_poly_evaluate_nmod_precomp(a->coeffs, a->length, c, c_precomp, a->mod);
            eval1_v2 = n_addmod(eval1_v2, 
                                _nmod_poly_evaluate_nmod_precomp(b->coeffs, b->length,
                                                               c, c_precomp, b->mod),
                                n);
#if FLINT_BITS == 64
            if (n <= UWORD(6148914691236517205))
#else // FLINT_BITS == 32
            if (n <= UWORD(1431655765))
#endif
            {
                eval1_v3 = _nmod_poly_evaluate_nmod_precomp_lazy(a->coeffs,
                                          a->length, c, c_precomp, a->mod);
                eval1_v3 = n_addmod(eval1_v3, 
                                    _nmod_poly_evaluate_nmod_precomp_lazy(b->coeffs,
                                                  b->length, c, c_precomp, b->mod),
                                    n);
            }
            else
                eval1_v3 = eval1_v1;
        }
        else
        {
            eval1_v2 = eval1_v1;
            eval1_v3 = eval1_v1;
        }

        nmod_poly_add(a, a, b);
        eval2_v1 = _nmod_poly_evaluate_nmod(a->coeffs, a->length, c, a->mod);
        if (a->mod.norm > 0)
        {
            ulong c_precomp = n_mulmod_precomp_shoup(c, n);
            eval2_v2 = _nmod_poly_evaluate_nmod_precomp(a->coeffs, a->length, c, c_precomp, a->mod);
#if FLINT_BITS == 64
            if (n <= UWORD(6148914691236517205))
#else // FLINT_BITS == 32
            if (n <= UWORD(1431655765))
#endif
                eval2_v3 = _nmod_poly_evaluate_nmod_precomp_lazy(a->coeffs, a->length, c, c_precomp, a->mod);
            else
                eval2_v3 = eval2_v1;
        }
        else
        {
            eval2_v2 = eval2_v1;
            eval2_v3 = eval2_v1;
        }


        result = (eval1 == eval2)
                 && (eval1 == eval1_v1) && (eval1 == eval1_v2) && (eval1 == eval1_v3)
                 && (eval2 == eval2_v1) && (eval2 == eval2_v2) && (eval2 == eval2_v3);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("eval1 = %wu, eval2 = %wu\n", eval1, eval2);
            flint_printf("eval1_v1 = %wu, eval2_v1 = %wu\n", eval1_v1, eval2_v1);
            flint_printf("eval1_v2 = %wu, eval2_v2 = %wu\n", eval1_v2, eval2_v2);
            flint_printf("eval1_v3 = %wu, eval2_v3 = %wu\n", eval1_v3, eval2_v3);
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(acopy);
        nmod_poly_clear(b);
    }

    TEST_FUNCTION_END(state);
}
