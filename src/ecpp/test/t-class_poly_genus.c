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

TEST_FUNCTION_START(ecpp_class_poly_genus, state)
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

    for (iter = 0; iter < 30 * flint_test_multiplier(); iter++)
    {
        ecpp_disc_struct * discs;
        slong ndiscs, k, g, o, tries;
        fmpz_t n, np, s;
        fmpz_mod_ctx_t ctx;
        fmpz_mod_poly_t F, H, R;
        fmpz_poly_t HZ;
        slong pstar[ECPP_DISC_MAXFAC + 1];
        fmpz sqrts[ECPP_DISC_MAXFAC + 1];
        const ecpp_disc_struct * d;

        ndiscs = ecpp_disc_table(&discs, primes, nprimes, 20000, 64, 64, 0.0);
        /* a discriminant with at least two genus characters */
        do
            k = n_randint(state, ndiscs);
        while (discs[k].g < 2);
        d = discs + k;

        fmpz_init(n); fmpz_init(np); fmpz_init(s);
        for (i = 0; i <= ECPP_DISC_MAXFAC; i++)
            fmpz_init(sqrts + i);

        /* a prime with all characters 1 */
        for (tries = 0; ; tries++)
        {
            int ok = 1;
            ulong n8;
            fmpz_randprime(n, state, 60 + n_randint(state, 100), 0);
            n8 = fmpz_fdiv_ui(n, 8);
            if (d->q0 == -4 && n8 % 4 != 1) ok = 0;
            if (d->q0 == 8 && !(n8 == 1 || n8 == 7)) ok = 0;
            if (d->q0 == -8 && !(n8 == 1 || n8 == 3)) ok = 0;
            for (i = 0; i < d->nfac && ok; i++)
                if (n_jacobi(fmpz_fdiv_ui(n, primes[d->fac[i]]), primes[d->fac[i]]) != 1)
                    ok = 0;
            if (ok)
                break;
        }

        fmpz_mod_ctx_init(ctx, n);
        g = 0;
        for (i = 0; i < d->nfac; i++)
        {
            ulong p = primes[d->fac[i]];
            pstar[g] = (p % 4 == 1) ? (slong) p : -(slong) p;
            fmpz_set_si(s, pstar[g]);
            fmpz_mod_set_fmpz(s, s, ctx);
            if (!fmpz_sqrtmod(sqrts + g, s, n))
                flint_abort();
            g++;
        }
        if (d->q0 != 1)
        {
            pstar[g] = d->q0;
            fmpz_set_si(s, d->q0 / 4);
            fmpz_mod_set_fmpz(s, s, ctx);
            if (!fmpz_sqrtmod(sqrts + g, s, n))
                flint_abort();
            fmpz_mod_add(sqrts + g, sqrts + g, sqrts + g, ctx);
            g++;
        }
        if (g != d->g)
        {
            flint_printf("FAIL: g mismatch for D = %wd\n", d->D);
            flint_abort();
        }

        fmpz_mod_poly_init(F, ctx);
        fmpz_mod_poly_init(H, ctx);
        fmpz_mod_poly_init(R, ctx);
        fmpz_poly_init(HZ);

        if (!ecpp_class_poly_genus(F, d->D, pstar, g, sqrts, ctx))
        {
            flint_printf("FAIL: genus factor not found for D = %wd (h = %wd, g = %wd)\n", d->D, d->h, d->g);
            flint_abort();
        }
        acb_modular_hilbert_class_poly(HZ, d->D);
        fmpz_mod_poly_set_fmpz_poly(H, HZ, ctx);
        o = d->h;
        for (i = 1; i < g; i++)
            o /= 2;
        if (fmpz_mod_poly_degree(F, ctx) != o)
        {
            flint_printf("FAIL: degree %wd, expected %wd (D = %wd)\n", fmpz_mod_poly_degree(F, ctx), o, d->D);
            flint_abort();
        }
        fmpz_mod_poly_rem(R, H, F, ctx);
        if (!fmpz_mod_poly_is_zero(R, ctx))
        {
            flint_printf("FAIL: factor does not divide H_D mod n (D = %wd, h = %wd, g = %wd)\n", d->D, d->h, d->g);
            flint_printf("n = "); fmpz_print(n); flint_printf("\n");
            flint_abort();
        }

        fmpz_mod_poly_clear(F, ctx);
        fmpz_mod_poly_clear(H, ctx);
        fmpz_mod_poly_clear(R, ctx);
        fmpz_poly_clear(HZ);
        fmpz_mod_ctx_clear(ctx);
        for (i = 0; i <= ECPP_DISC_MAXFAC; i++)
            fmpz_clear(sqrts + i);
        fmpz_clear(n); fmpz_clear(np); fmpz_clear(s);
        flint_free(discs);
    }

    TEST_FUNCTION_END(state);
}
