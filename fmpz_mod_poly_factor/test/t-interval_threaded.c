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

#undef ulong

#include <gmp.h>
#include <pthread.h>

#define ulong mp_limb_t

#include "flint.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    FLINT_TEST_INIT(state);
    
    flint_printf("interval_threaded....");
    fflush(stdout);

#if HAVE_PTHREAD && (HAVE_TLS || FLINT_REENTRANT)

    /* no aliasing */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        fmpz_mod_poly_t a, b, c, cinv, d, *e, * tmp;
        fmpz_t p;
        slong j, num_threads, l;
        fmpz_mod_poly_interval_poly_arg_t * args1;
        pthread_t *threads;

        flint_set_num_threads(1 + n_randint(state, 3));

        num_threads = flint_get_num_threads();

        l = n_randint(state, 20) + 1;
        threads = flint_malloc(sizeof(pthread_t) * num_threads);
        e = flint_malloc(sizeof(fmpz_mod_poly_struct) * num_threads);
        tmp = flint_malloc(sizeof(fmpz_mod_poly_struct) * l);
        args1 = flint_malloc(num_threads *
                             sizeof(fmpz_mod_poly_interval_poly_arg_t));

        fmpz_init(p);
        fmpz_set_ui(p, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(a, p);
        fmpz_mod_poly_init(b, p);
        fmpz_mod_poly_init(c, p);
        fmpz_mod_poly_init(cinv, p);
        fmpz_mod_poly_init(d, p);
        for (j = 0; j < l; j++)
            fmpz_mod_poly_init(tmp[j], p);
        for (j = 0; j < num_threads; j++)
            fmpz_mod_poly_init(e[j], p);

        fmpz_mod_poly_randtest_not_zero(a, state, n_randint(state, 20) + 1);
        do
        {
            fmpz_mod_poly_randtest_not_zero(c, state, n_randint(state, 20) + 1);
        } while (c->length < 3);

        fmpz_mod_poly_rem(a, a, c);

        for (j = 0; j < l; j++)
            fmpz_mod_poly_randtest_not_zero(tmp[j], state, n_randint(state, 20) + 1);

        fmpz_mod_poly_reverse(cinv, c, c->length);
        fmpz_mod_poly_inv_series_newton(cinv, cinv, c->length);

        fmpz_mod_poly_set_ui(b, UWORD(1));
        for (j = l - 1; j >= 0; j--)
        {
            fmpz_mod_poly_rem(d, tmp[j], c);
            fmpz_mod_poly_sub(d, a, d);
            fmpz_mod_poly_mulmod_preinv(b, d, b, c, cinv);
        }


        for (j = 0; j < num_threads; j++)
        {
            fmpz_mod_poly_fit_length(e[j], c->length - 1);
            _fmpz_mod_poly_set_length(e[j], c->length - 1);
            _fmpz_vec_zero(e[j]->coeffs, c->length - 1);
            args1[j].baby = *tmp;
            args1[j].res = *e[j];
            args1[j].H = *a;
            args1[j].v = *c;
            args1[j].vinv = *cinv;
            args1[j].m = l;

            pthread_create(&threads[j], NULL, _fmpz_mod_poly_interval_poly_worker, &args1[j]);
        }
        for (j = 0; j < num_threads; j++)
            pthread_join(threads[j], NULL);
        for (j = 0; j < num_threads; j++)
            _fmpz_mod_poly_normalise(e[j]);

        for (j = 0; j < num_threads; j++)
        {
            if (!fmpz_mod_poly_equal(b, e[j]))
            {
                flint_printf("j: %wd\n", j);
                flint_printf("FAIL (interval_poly):\n");
                flint_printf("b:\n"); fmpz_mod_poly_print(b); flint_printf("\n");
                flint_printf("c:\n"); fmpz_mod_poly_print(c); flint_printf("\n");
                flint_printf("e[j]:\n"); fmpz_mod_poly_print(e[j]); flint_printf("\n");
                abort();
            }
        }

        fmpz_clear(p);
        fmpz_mod_poly_clear(a);
        fmpz_mod_poly_clear(b);
        fmpz_mod_poly_clear(c);
        fmpz_mod_poly_clear(cinv);
        fmpz_mod_poly_clear(d);
        for (j = 0; j < num_threads; j++)
            fmpz_mod_poly_clear(e[j]);
        for (j = 0; j < l; j++)
            fmpz_mod_poly_clear(tmp[j]);
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

