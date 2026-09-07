/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.1

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "acb_modular.h"
#include "ecpp.h"

TEST_FUNCTION_START(ecpp_class_poly_tower, state)
{
    slong iter;
    ulong primes[64];
    slong nprimes = 0, i;
    n_primes_t it;

    n_primes_init(it);
    n_primes_next(it);
    for (i = 0; i < 64; i++)
        primes[nprimes++] = n_primes_next(it);
    n_primes_clear(it);

    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        ecpp_disc_struct * discs;
        slong ndiscs, k, tries;
        fmpz_t n, j, t, v, sq;
        fmpz_mod_ctx_t ctx;
        fmpz_poly_t HZ;
        fmpz_mod_poly_t H;
        const ecpp_disc_struct * d;

        ndiscs = ecpp_disc_table(&discs, primes, nprimes, 40000, 48, 48, 0.0);
        do
            k = n_randint(state, ndiscs);
        while (discs[k].h < 2);
        d = discs + k;

        fmpz_init(n); fmpz_init(j); fmpz_init(t); fmpz_init(v); fmpz_init(sq);

        /* a prime represented by the principal form: 4n = t^2 + |D| v^2 */
        for (tries = 0; ; tries++)
        {
            fmpz_randprime(n, state, 60 + n_randint(state, 100), 0);
            fmpz_set_si(sq, d->D);
            fmpz_mod(sq, sq, n);
            if (fmpz_jacobi(sq, n) != 1)
                continue;
            fmpz_sqrtmod(sq, sq, n);
            if (ecpp_cornacchia(t, v, n, d->D, sq))
                break;
        }

        fmpz_mod_ctx_init(ctx, n);
        fmpz_poly_init(HZ);
        fmpz_mod_poly_init(H, ctx);

        if (!ecpp_class_poly_tower(j, d->D, state, ctx))
        {
            flint_printf("FAIL: no root found for D = %wd (h = %wd)\n", d->D, d->h);
            flint_abort();
        }
        acb_modular_hilbert_class_poly(HZ, d->D);
        fmpz_mod_poly_set_fmpz_poly(H, HZ, ctx);
        fmpz_mod_poly_evaluate_fmpz(t, H, j, ctx);
        if (!fmpz_is_zero(t))
        {
            flint_printf("FAIL: not a root of H_D for D = %wd (h = %wd)\n", d->D, d->h);
            flint_printf("n = "); fmpz_print(n); flint_printf("\n");
            flint_abort();
        }

        /* odd D with v even: the order of conductor 2 with Weber's f */
        if (d->D % 2 != 0 && fmpz_is_even(v))
        {
            if (!ecpp_class_poly_tower(j, 4 * d->D, state, ctx))
            {
                flint_printf("FAIL: no root found for 4D, D = %wd (h = %wd)\n", d->D, d->h);
                flint_abort();
            }
            acb_modular_hilbert_class_poly(HZ, 4 * d->D);
            fmpz_mod_poly_set_fmpz_poly(H, HZ, ctx);
            fmpz_mod_poly_evaluate_fmpz(t, H, j, ctx);
            if (!fmpz_is_zero(t))
            {
                flint_printf("FAIL: not a root of H_{4D} for D = %wd (h = %wd)\n", d->D, d->h);
                flint_abort();
            }
        }

        fmpz_mod_poly_clear(H, ctx);
        fmpz_poly_clear(HZ);
        fmpz_mod_ctx_clear(ctx);
        fmpz_clear(n); fmpz_clear(j); fmpz_clear(t); fmpz_clear(v); fmpz_clear(sq);
        flint_free(discs);
    }

    TEST_FUNCTION_END(state);
}
