/*
    Copyright (C) 2014 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#undef ulong
#define ulong ulongxx/* interferes with system includes */

#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

#undef ulong

#include <gmp.h>

#define ulong mp_limb_t

#include "flint.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int iter;
    FLINT_TEST_INIT(state);

    flint_printf("interval_threaded....");
    fflush(stdout);

#if HAVE_PTHREAD && (HAVE_TLS || FLINT_REENTRANT)

    for (iter = 0; iter < 20*flint_test_multiplier(); iter++)
    {
        nmod_poly_t a, b, c, cinv, d, *e, * tmp;
        mp_limb_t modulus;
        slong j, num_threads, l;
        nmod_poly_interval_poly_arg_t * args1;
        pthread_t *threads;

        flint_set_num_threads(1 + n_randint(state, 3));

        num_threads = flint_get_num_threads();

        l = n_randint(state, 20) + 1;
        threads = flint_malloc(sizeof(pthread_t) * num_threads);
        e = flint_malloc(sizeof(nmod_poly_struct) * num_threads);
        tmp = flint_malloc(sizeof(nmod_poly_struct) * l);
        args1 = flint_malloc(num_threads *
                             sizeof(nmod_poly_interval_poly_arg_t));

        modulus = n_randtest_prime(state, 0);

        nmod_poly_init(a, modulus);
        nmod_poly_init(b, modulus);
        nmod_poly_init(c, modulus);
        nmod_poly_init(cinv, modulus);
        nmod_poly_init(d, modulus);
        for (j = 0; j < l; j++)
            nmod_poly_init(tmp[j], modulus);
        for (j = 0; j < num_threads; j++)
            nmod_poly_init(e[j], modulus);

        nmod_poly_randtest_not_zero(a, state, n_randint(state, 20) + 1);
        do
        {
            nmod_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1);
        } while (c->length < 3);

        nmod_poly_rem(a, a, c);

        for (j = 0; j < l; j++)
            nmod_poly_randtest_not_zero(tmp[j], state, n_randint(state, 20) + 1);

        nmod_poly_reverse(cinv, c, c->length);
        nmod_poly_inv_series(cinv, cinv, c->length);

        nmod_poly_one(b);
        for (j = l - 1; j >= 0; j--)
        {
            nmod_poly_rem(d, tmp[j], c);
            nmod_poly_sub(d, a, d);
            nmod_poly_mulmod_preinv(b, d, b, c, cinv);
        }


        for (j = 0; j < num_threads; j++)
        {
            nmod_poly_fit_length(e[j], c->length - 1);
            _nmod_poly_set_length(e[j], c->length - 1);
            _nmod_vec_zero(e[j]->coeffs, c->length - 1);
            args1[j].baby = *tmp;
            args1[j].res = *e[j];
            args1[j].H = *a;
            args1[j].v = *c;
            args1[j].vinv = *cinv;
            args1[j].m = l;

            pthread_create(&threads[j], NULL, _nmod_poly_interval_poly_worker, &args1[j]);
        }
        for (j = 0; j < num_threads; j++)
            pthread_join(threads[j], NULL);
        for (j = 0; j < num_threads; j++)
            _nmod_poly_normalise(e[j]);

        for (j = 0; j < num_threads; j++)
        {
            if (!nmod_poly_equal(b, e[j]))
            {
                flint_printf("j: %wd\n", j);
                flint_printf("FAIL (interval_poly):\n");
                flint_printf("b:\n"); nmod_poly_print(b); flint_printf("\n");
                flint_printf("c:\n"); nmod_poly_print(c); flint_printf("\n");
                flint_printf("e[j]:\n"); nmod_poly_print(e[j]); flint_printf("\n");
                abort();
            }
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(cinv);
        nmod_poly_clear(d);
        for (j = 0; j < num_threads; j++)
            nmod_poly_clear(e[j]);
        for (j = 0; j < l; j++)
            nmod_poly_clear(tmp[j]);
        flint_free(e);
        flint_free(tmp);
        flint_free(args1);
        flint_free(threads);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;

#else

   FLINT_TEST_CLEANUP(state);

   flint_printf("SKIPPED\n");
   return 0;

#endif

}

