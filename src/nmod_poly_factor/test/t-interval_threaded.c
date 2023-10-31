/*
    Copyright (C) 2014 Martin Lee
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "thread_support.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_poly_factor.h"

TEST_FUNCTION_START(nmod_poly_factor_interval_threaded, state)
{
#if FLINT_USES_PTHREAD && (FLINT_USES_TLS || FLINT_REENTRANT)
    int iter;

    for (iter = 0; iter < 20*flint_test_multiplier(); iter++)
    {
        nmod_poly_t a, b, c, cinv, d;
	nmod_poly_struct * tmp;
	nmod_poly_struct * e;
        mp_limb_t modulus;
        slong j, num_threads, l;
        nmod_poly_interval_poly_arg_t * args1;
        thread_pool_handle * threads;

        flint_set_num_threads(1 + n_randint(state, 3));

        num_threads = flint_request_threads(&threads,
			                           FLINT_DEFAULT_THREAD_LIMIT);

        l = n_randint(state, 20) + 1;

        e = (nmod_poly_struct *)
		flint_malloc(sizeof(nmod_poly_struct)*(num_threads + 1));
        tmp = (nmod_poly_struct *)
		flint_malloc(sizeof(nmod_poly_struct)*l);
        args1 = (nmod_poly_interval_poly_arg_t *)
		flint_malloc((num_threads + 1)*
                             sizeof(nmod_poly_interval_poly_arg_t));

        modulus = n_randtest_prime(state, 0);

        nmod_poly_init(a, modulus);
        nmod_poly_init(b, modulus);
        nmod_poly_init(c, modulus);
        nmod_poly_init(cinv, modulus);
        nmod_poly_init(d, modulus);

	for (j = 0; j < l; j++)
            nmod_poly_init(tmp + j, modulus);

	for (j = 0; j < num_threads + 1; j++)
            nmod_poly_init(e + j, modulus);

        nmod_poly_randtest_not_zero(a, state, n_randint(state, 20) + 1);

	do
        {
            nmod_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1);
        } while (c->length < 3);

        nmod_poly_rem(a, a, c);

        for (j = 0; j < l; j++)
            nmod_poly_randtest_not_zero(tmp + j, state, n_randint(state, 20) + 1);

        nmod_poly_reverse(cinv, c, c->length);
        nmod_poly_inv_series(cinv, cinv, c->length);

        nmod_poly_one(b);
	for (j = l - 1; j >= 0; j--)
        {
            nmod_poly_rem(d, tmp + j, c);
            nmod_poly_sub(d, a, d);
            nmod_poly_mulmod_preinv(b, d, b, c, cinv);
        }

        for (j = 0; j < num_threads + 1; j++)
        {
            nmod_poly_fit_length(e + j, c->length - 1);
            _nmod_poly_set_length(e + j, c->length - 1);
            _nmod_vec_zero(e[j].coeffs, c->length - 1);

	    args1[j].baby = tmp;
            args1[j].res = e + j;
            args1[j].H = a;
            args1[j].v = c;
            args1[j].vinv = cinv;
            args1[j].m = l;
	    args1[j].tmp = _nmod_vec_init(c->length - 1);
        }

	for (j = 0; j < num_threads; j++)
		thread_pool_wake(global_thread_pool, threads[j], 0,
				   _nmod_poly_interval_poly_worker, &args1[j]);

	_nmod_poly_interval_poly_worker(&args1[num_threads]);

        for (j = 0; j < num_threads; j++)
            thread_pool_wait(global_thread_pool, threads[j]);

	for (j = 0; j < num_threads + 1; j++)
        {
            _nmod_poly_normalise(e + j);
	    _nmod_vec_clear(args1[j].tmp);
	}

        for (j = 0; j < num_threads + 1; j++)
        {
            if (!nmod_poly_equal(b, e + j))
            {
                flint_printf("j: %wd\n", j);
                flint_printf("FAIL (interval_poly):\n");
                flint_printf("b:\n"); nmod_poly_print(b); flint_printf("\n");
                flint_printf("c:\n"); nmod_poly_print(c); flint_printf("\n");
                flint_printf("e[j]:\n"); nmod_poly_print(e + j); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

	flint_give_back_threads(threads, num_threads);

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(cinv);
        nmod_poly_clear(d);

        for (j = 0; j < num_threads + 1; j++)
            nmod_poly_clear(e + j);

        for (j = 0; j < l; j++)
            nmod_poly_clear(tmp + j);

	flint_free(e);
        flint_free(tmp);
        flint_free(args1);
    }

    TEST_FUNCTION_END(state);
#else
    TEST_FUNCTION_END_SKIPPED(state);
#endif
}
