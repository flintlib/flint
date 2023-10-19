/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Martin Lee
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
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

TEST_FUNCTION_START(fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv_threaded, state)
{
#if FLINT_USES_PTHREAD && (FLINT_USES_TLS || FLINT_REENTRANT)
    int i;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* check precomputation */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, c, cinv;
        fmpz_mod_poly_struct * tmp;
        fmpz_t p;
        fmpz_mat_t B;
        fmpz_mat_struct * C;
        slong j, num_threads;
        fmpz_mod_poly_matrix_precompute_arg_t * args1;
        thread_pool_handle * threads;

        flint_set_num_threads(1 + n_randint(state, 3));

        num_threads = flint_request_threads(&threads,
			                           FLINT_DEFAULT_THREAD_LIMIT);

        tmp = (fmpz_mod_poly_struct *)
                  flint_malloc(sizeof(fmpz_mod_poly_struct)*(num_threads + 1));

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_init(cinv, ctx);

        for (j = 0; j < num_threads + 1; j++)
            fmpz_mod_poly_init(tmp + j, ctx);

        fmpz_mod_poly_randtest_not_zero(a, state, n_randint(state, 20) + 1, ctx);
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 20) + 1, ctx);
        do
        {
            fmpz_mod_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1, ctx);
        } while (c->length < 2);

        fmpz_mod_poly_reverse(cinv, c, c->length, ctx);
        fmpz_mod_poly_inv_series(cinv, cinv, c->length, ctx);

        fmpz_mat_init(B, n_sqrt(c->length - 1) + 1, c->length - 1);
        fmpz_mod_poly_precompute_matrix(B, b, c, cinv, ctx);

        args1 = (fmpz_mod_poly_matrix_precompute_arg_t *)
               flint_malloc(sizeof(fmpz_mod_poly_matrix_precompute_arg_t)
                                                           *(num_threads + 1));
        C = (fmpz_mat_struct* )
                       flint_malloc(sizeof(fmpz_mat_struct)*(num_threads + 1));

        for (j = 0; j < num_threads + 1; j++)
        {
            fmpz_mat_init(C + j, n_sqrt(c->length - 1) + 1, c->length - 1);

    	    fmpz_mod_poly_set(tmp + j, b, ctx);
            fmpz_mod_poly_rem(tmp + j, tmp + j, c, ctx);

            if (tmp[j].length < c->length - 1)
            {
                fmpz_mod_poly_fit_length(tmp + j, c->length - 1, ctx);
                _fmpz_vec_zero(tmp[j].coeffs + tmp[j].length,
                                                    c->length - 1 - b->length);
            }

            args1[j].A        = C + j;
            args1[j].poly1    = tmp + j;
            args1[j].poly2    = c;
            args1[j].poly2inv = cinv;
            args1[j].ctx      = ctx;
        }

        for (j = 1; j < num_threads + 1; j++)
                thread_pool_wake(global_thread_pool, threads[j - 1], 0,
                           _fmpz_mod_poly_precompute_matrix_worker, &args1[j]);

        _fmpz_mod_poly_precompute_matrix_worker(&args1[0]);

        for (j = 0; j < num_threads; j++)
        {
            thread_pool_wait(global_thread_pool, threads[j]);
        }

        for (j = 0; j < num_threads + 1; j++)
        {
            if (!fmpz_mat_equal(B, C + j))
            {
                flint_printf("FAIL (precomputation):\n");
                flint_printf("B:\n"); fmpz_mat_print(B); flint_printf("\n");
                flint_printf("C[j]:\n"); fmpz_mat_print(C + j); flint_printf("\n");
                flint_printf("a:\n"); fmpz_mod_poly_print(a, ctx); flint_printf("\n");
                flint_printf("b:\n"); fmpz_mod_poly_print(b, ctx); flint_printf("\n");
                flint_printf("c:\n"); fmpz_mod_poly_print(c, ctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        flint_give_back_threads(threads, num_threads);

        fmpz_clear(p);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mat_clear(B);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_mod_poly_clear(cinv, ctx);

        for (j = 0; j < num_threads + 1; j++)
        {
            fmpz_mod_poly_clear(tmp + j, ctx);
            fmpz_mat_clear(C + j);
        }

        flint_free(C);
        flint_free(tmp);
        flint_free(args1);
    }

    /* check composition */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, c, cinv, d;
        fmpz_mod_poly_struct * res;
        fmpz_t p;
        fmpz_mat_t B;
        slong j, num_threads;
        fmpz_mod_poly_compose_mod_precomp_preinv_arg_t * args1;
        thread_pool_handle * threads;

        flint_set_num_threads(1 + n_randint(state, 3));

        num_threads = flint_request_threads(&threads,
			                           FLINT_DEFAULT_THREAD_LIMIT);

        res = (fmpz_mod_poly_struct* )
            flint_malloc(sizeof(fmpz_mod_poly_struct)*(num_threads + 1));

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_init(cinv, ctx);
        fmpz_mod_poly_init(d, ctx);

        for (j = 0; j < num_threads + 1; j++)
            fmpz_mod_poly_init(res + j, ctx);

        fmpz_mod_poly_randtest_not_zero(a, state, n_randint(state, 20) + 1, ctx);
        fmpz_mod_poly_randtest_not_zero(b, state, n_randint(state, 20) + 1, ctx);
        do
        {
            fmpz_mod_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1, ctx);
        } while (c->length < 2);

        fmpz_mod_poly_reverse(cinv, c, c->length, ctx);
        fmpz_mod_poly_inv_series(cinv, cinv, c->length, ctx);

        fmpz_mat_init(B, n_sqrt(c->length - 1) + 1, c->length - 1);
        fmpz_mod_poly_precompute_matrix(B, b, c, cinv, ctx);

        fmpz_mod_poly_rem(a, a, c, ctx);
        fmpz_mod_poly_compose_mod(d, a, b, c, ctx);

        args1 = (fmpz_mod_poly_compose_mod_precomp_preinv_arg_t *)
                   flint_malloc((num_threads + 1)*
                        sizeof(fmpz_mod_poly_compose_mod_precomp_preinv_arg_t));

        for (j = 0; j < num_threads + 1; j++)
        {
            fmpz_mod_poly_fit_length(res + j, c->length - 1, ctx);
            _fmpz_mod_poly_set_length(res + j, c->length - 1);

            args1[j].A        = B;
            args1[j].res      = res + j;
            args1[j].poly1    = a;
            args1[j].poly3    = c;
            args1[j].poly3inv = cinv;
            args1[j].ctx      = ctx;
        }

        for (j = 1; j < num_threads + 1; j++)
          thread_pool_wake(global_thread_pool, threads[j - 1], 0,
              _fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv_worker,
		                                                    &args1[j]);

        _fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv_worker(&args1[0]);
        _fmpz_mod_poly_normalise(res + 0);

        for (j = 1; j < num_threads + 1; j++)
        {
            thread_pool_wait(global_thread_pool, threads[j - 1]);
            _fmpz_mod_poly_normalise(res + j);
        }

        for (j = 0; j < num_threads + 1; j++)
        {
            if (!fmpz_mod_poly_equal(d, res + j, ctx))
            {
                flint_printf("FAIL (composition):\n");
                flint_printf("res[j]:\n"); fmpz_mod_poly_print(res + j, ctx); flint_printf("\n");
                flint_printf("d:\n"); fmpz_mod_poly_print(d, ctx); flint_printf("\n");
                flint_printf("a:\n"); fmpz_mod_poly_print(a, ctx); flint_printf("\n");
                flint_printf("b:\n"); fmpz_mod_poly_print(b, ctx); flint_printf("\n");
                flint_printf("c:\n"); fmpz_mod_poly_print(c, ctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        flint_give_back_threads(threads, num_threads);

        fmpz_clear(p);
        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mat_clear(B);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_mod_poly_clear(cinv, ctx);
        fmpz_mod_poly_clear(d, ctx);

        for (j = 0; j < num_threads + 1; j++)
            fmpz_mod_poly_clear(res + j, ctx);

        flint_free(res);
        flint_free(args1);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
#else
    TEST_FUNCTION_END_SKIPPED(state);
#endif
}
