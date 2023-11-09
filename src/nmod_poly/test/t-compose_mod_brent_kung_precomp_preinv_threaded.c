/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013, 2014 Martin Lee
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "thread_support.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_poly.h"

TEST_FUNCTION_START(nmod_poly_compose_mod_brent_kung_precomp_preinv_threaded, state)
{
#if FLINT_USES_PTHREAD && (FLINT_USES_TLS || FLINT_REENTRANT)
    int i;

    /* check precomputation */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, cinv;
	nmod_poly_struct * tmp;
        nmod_mat_t B;
	nmod_mat_struct * C;
        mp_limb_t m = n_randtest_prime(state, 0);
        slong j, num_threads;
        nmod_poly_matrix_precompute_arg_t * args1;
        thread_pool_handle * threads;

        flint_set_num_threads(1 + n_randint(state, 3));

        num_threads = flint_request_threads(&threads, FLINT_DEFAULT_THREAD_LIMIT);

        tmp = (nmod_poly_struct *) flint_malloc(sizeof(nmod_poly_struct)*(num_threads + 1));

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        nmod_poly_init(cinv, m);

        for (j = 0; j < num_threads + 1; j++)
            nmod_poly_init(tmp + j, m);

        nmod_poly_randtest_not_zero(a, state, 1 + n_randint(state, 20));
        nmod_poly_randtest_not_zero(b, state, 1 + n_randint(state, 20));
        do
        {
            nmod_poly_randtest_not_zero(c, state, 1 + n_randint(state, 20));
        } while (c->length < 2);

        nmod_poly_reverse(cinv, c, c->length);
        nmod_poly_inv_series(cinv, cinv, c->length);

        nmod_mat_init(B, n_sqrt (c->length - 1) + 1, c->length - 1, m);
        nmod_poly_precompute_matrix(B, b, c, cinv);

        args1 = (nmod_poly_matrix_precompute_arg_t *)
		 flint_malloc(sizeof(nmod_poly_matrix_precompute_arg_t)
                             *(num_threads + 1));
        C = (nmod_mat_struct *)
	               flint_malloc(sizeof(nmod_mat_struct)*(num_threads + 1));

        for (j = 0; j < num_threads + 1; j++)
        {
            nmod_mat_init(C + j, n_sqrt(c->length - 1) + 1, c->length - 1, m);
            nmod_poly_set(tmp + j, b);
            nmod_poly_rem(tmp + j, tmp + j, c);

	    if (tmp[j].length < c->length - 1)
            {
                nmod_poly_fit_length(tmp + j, c->length - 1);
                _nmod_vec_zero(tmp[j].coeffs + tmp[j].length,
                               c->length - 1 - b->length);
            }

            args1[j].A        = C + j;
            args1[j].poly1    = tmp + j;
            args1[j].poly2    = c;
            args1[j].poly2inv = cinv;
        }

	for (j = 1; j < num_threads + 1; j++)
        {
	    thread_pool_wake(global_thread_pool, threads[j - 1], 0,
                           _nmod_poly_precompute_matrix_worker, &args1[j]);
        }

	_nmod_poly_precompute_matrix_worker(&args1[0]);

	for (j = 1; j < num_threads + 1; j++)
            thread_pool_wait(global_thread_pool, threads[j - 1]);

        for (j = 0; j < num_threads + 1; j++)
        {
            if (!nmod_mat_equal(B, C + j))
            {
                flint_printf("FAIL (precomputation):\n");
                flint_printf("B:\n"); nmod_mat_print_pretty(B); flint_printf("\n");
                flint_printf("C[j]:\n"); nmod_mat_print_pretty(C + j); flint_printf("\n");
                flint_printf("a:\n"); nmod_poly_print(a); flint_printf("\n");
                flint_printf("b:\n"); nmod_poly_print(b); flint_printf("\n");
                flint_printf("c:\n"); nmod_poly_print(c); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

	flint_give_back_threads(threads, num_threads);

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_mat_clear (B);
        nmod_poly_clear(c);
        nmod_poly_clear(cinv);

	for (j = 0; j < num_threads + 1; j++)
        {
            nmod_poly_clear(tmp + j);
            nmod_mat_clear(C + j);
        }

	flint_free(C);
        flint_free(tmp);
        flint_free(args1);
    }

    /* check composition */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, cinv, d;
	nmod_poly_struct * res;
        nmod_mat_t B;
        mp_limb_t m = n_randtest_prime(state, 0);
        slong j, num_threads;
        nmod_poly_compose_mod_precomp_preinv_arg_t * args1;
        thread_pool_handle * threads;

        flint_set_num_threads(1 + n_randint(state, 3));

        num_threads = flint_request_threads(&threads, FLINT_DEFAULT_THREAD_LIMIT);

        res = (nmod_poly_struct *)
		flint_malloc(sizeof(nmod_poly_struct)*(num_threads + 1));

        nmod_poly_init(a, m);
        nmod_poly_init(b, m);
        nmod_poly_init(c, m);
        nmod_poly_init(cinv, m);
        nmod_poly_init(d, m);

	for (j = 0; j < num_threads + 1; j++)
            nmod_poly_init(res + j, m);

        nmod_poly_randtest_not_zero(a, state, 1 + n_randint(state, 20));
        nmod_poly_randtest_not_zero(b, state, 1 + n_randint(state, 20));
        do
        {
            nmod_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1);
        } while (c->length < 2);

        nmod_poly_reverse(cinv, c, c->length);
        nmod_poly_inv_series(cinv, cinv, c->length);

        nmod_mat_init(B, n_sqrt(c->length - 1) + 1, c->length - 1, m);
        nmod_poly_precompute_matrix(B, b, c, cinv);

        nmod_poly_rem(a, a, c);
        nmod_poly_compose_mod(d, a, b, c);

        args1 = (nmod_poly_compose_mod_precomp_preinv_arg_t *)
		flint_malloc((num_threads + 1)*
                        sizeof(nmod_poly_compose_mod_precomp_preinv_arg_t));

        for (j = 0; j < num_threads + 1; j++)
        {
            nmod_poly_fit_length(res + j, c->length - 1);
            _nmod_poly_set_length(res + j, c->length - 1);
            flint_mpn_zero(res[j].coeffs, c->length - 1);

	    args1[j].A        = B;
            args1[j].res      = res + j;
            args1[j].poly1    = a;
            args1[j].poly3    = c;
            args1[j].poly3inv = cinv;
	}

        for (j = 1; j < num_threads + 1; j++)
	{
            thread_pool_wake(global_thread_pool, threads[j - 1], 0,
                _nmod_poly_compose_mod_brent_kung_precomp_preinv_worker, &args1[j]);
        }

	_nmod_poly_compose_mod_brent_kung_precomp_preinv_worker(&args1[0]);

	_nmod_poly_normalise(res + 0);

	for (j = 1; j < num_threads + 1; j++)
        {
            thread_pool_wait(global_thread_pool, threads[j - 1]);
            _nmod_poly_normalise(res + j);
        }

        for (j = 0; j < num_threads + 1; j++)
        {
            if (!nmod_poly_equal(d, res + j))
            {
                flint_printf("FAIL (composition):\n");
                flint_printf("res[j]:\n"); nmod_poly_print(res + j); flint_printf("\n");
                flint_printf("d:\n"); nmod_poly_print(d); flint_printf("\n");
                flint_printf("a:\n"); nmod_poly_print(a); flint_printf("\n");
                flint_printf("b:\n"); nmod_poly_print(b); flint_printf("\n");
                flint_printf("c:\n"); nmod_poly_print(c); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

	flint_give_back_threads(threads, num_threads);

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_mat_clear (B);
        nmod_poly_clear(c);
        nmod_poly_clear(cinv);
        nmod_poly_clear(d);

        for (j = 0; j < num_threads + 1; j++)
            nmod_poly_clear(res + j);

	flint_free(res);
        flint_free(args1);
    }

    TEST_FUNCTION_END(state);

#else
    TEST_FUNCTION_END_SKIPPED(state);
#endif
}
