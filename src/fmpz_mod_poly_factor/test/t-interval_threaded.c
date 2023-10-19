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
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mod_poly_factor.h"

TEST_FUNCTION_START(fmpz_mod_poly_factor_interval_threaded, state)
{
#if FLINT_USES_PTHREAD && (FLINT_USES_TLS || FLINT_REENTRANT)
    int i;
    fmpz_mod_ctx_t ctx;

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* no aliasing */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, c, cinv, d;
        fmpz_mod_poly_struct * tmp;
        fmpz_mod_poly_struct * e;
        fmpz_t p;
        slong j, num_threads, l;
        fmpz_mod_poly_interval_poly_arg_t * args1;
        thread_pool_handle * threads;

        flint_set_num_threads(1 + n_randint(state, 3));

        num_threads = flint_request_threads(&threads,
			                           FLINT_DEFAULT_THREAD_LIMIT);

        l = n_randint(state, 20) + 1;

        e = (fmpz_mod_poly_struct *)
             flint_malloc(sizeof(fmpz_mod_poly_struct)*(num_threads + 1));
        tmp = (fmpz_mod_poly_struct *)
             flint_malloc(sizeof(fmpz_mod_poly_struct)*l);
        args1 = (fmpz_mod_poly_interval_poly_arg_t *)
             flint_malloc((num_threads + 1)*
                             sizeof(fmpz_mod_poly_interval_poly_arg_t));

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));
        fmpz_mod_ctx_set_modulus(ctx, p);

        fmpz_mod_poly_init(a, ctx);
        fmpz_mod_poly_init(b, ctx);
        fmpz_mod_poly_init(c, ctx);
        fmpz_mod_poly_init(cinv, ctx);
        fmpz_mod_poly_init(d, ctx);

        for (j = 0; j < l; j++)
            fmpz_mod_poly_init(tmp + j, ctx);

        for (j = 0; j < num_threads + 1; j++)
            fmpz_mod_poly_init(e + j, ctx);

        fmpz_mod_poly_randtest_not_zero(a, state, n_randint(state, 20) + 1, ctx);

        do
        {
            fmpz_mod_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1, ctx);
        } while (c->length < 3);

        fmpz_mod_poly_rem(a, a, c, ctx);

        for (j = 0; j < l; j++)
            fmpz_mod_poly_randtest_not_zero(tmp + j, state, n_randint(state, 20) + 1, ctx);

        fmpz_mod_poly_reverse(cinv, c, c->length, ctx);
        fmpz_mod_poly_inv_series(cinv, cinv, c->length, ctx);

        fmpz_mod_poly_set_ui(b, 1, ctx);

        for (j = l - 1; j >= 0; j--)
        {
            fmpz_mod_poly_rem(d, tmp + j, c, ctx);
            fmpz_mod_poly_sub(d, a, d, ctx);
            fmpz_mod_poly_mulmod_preinv(b, d, b, c, cinv, ctx);
        }

        for (j = 0; j < num_threads + 1; j++)
        {
            fmpz_mod_poly_fit_length(e + j, c->length - 1, ctx);
            _fmpz_mod_poly_set_length(e + j, c->length - 1);
            _fmpz_vec_zero(e[j].coeffs, c->length - 1);

            args1[j].baby = tmp;
            args1[j].res = e + j;
            args1[j].H = a;
            args1[j].v = c;
            args1[j].vinv = cinv;
            args1[j].m = l;
            args1[j].tmp = _fmpz_vec_init(c->length - 1);
            args1[j].ctx = ctx;
        }

        for (j = 0; j < num_threads; j++)
            thread_pool_wake(global_thread_pool, threads[j], 0,
				   _fmpz_mod_poly_interval_poly_worker, &args1[j]);

        _fmpz_mod_poly_interval_poly_worker(&args1[num_threads]);

        for (j = 0; j < num_threads; j++)
            thread_pool_wait(global_thread_pool, threads[j]);

        for (j = 0; j < num_threads + 1; j++)
        {
            _fmpz_mod_poly_normalise(e + j);
            _fmpz_vec_clear(args1[j].tmp, c->length - 1);
        }

        for (j = 0; j < num_threads + 1; j++)
        {
            if (!fmpz_mod_poly_equal(b, e + j, ctx))
            {
                flint_printf("j: %wd\n", j);
                flint_printf("FAIL (interval_poly):\n");
                flint_printf("b:\n"); fmpz_mod_poly_print(b, ctx); flint_printf("\n");
                flint_printf("c:\n"); fmpz_mod_poly_print(c, ctx); flint_printf("\n");
                flint_printf("e[j]:\n"); fmpz_mod_poly_print(e + j, ctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        flint_give_back_threads(threads, num_threads);

        fmpz_clear(p);

        fmpz_mod_poly_clear(a, ctx);
        fmpz_mod_poly_clear(b, ctx);
        fmpz_mod_poly_clear(c, ctx);
        fmpz_mod_poly_clear(cinv, ctx);
        fmpz_mod_poly_clear(d, ctx);

        for (j = 0; j < num_threads + 1; j++)
            fmpz_mod_poly_clear(e + j, ctx);

        for (j = 0; j < l; j++)
            fmpz_mod_poly_clear(tmp + j, ctx);

        flint_free(e);
        flint_free(tmp);
        flint_free(args1);
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
#else
    TEST_FUNCTION_END_SKIPPED(state);
#endif
}
