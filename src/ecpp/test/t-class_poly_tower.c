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

/* a prime of the given size with n = residue mod modulus that splits
   completely in the Hilbert class field of D (4n = t^2 + |D| v^2), with
   v of the requested parity (or any, vpar < 0) */
static void
_split_prime(fmpz_t n, fmpz_t v, slong D, slong bits, ulong modulus, ulong residue,
                                                int vpar, flint_rand_t state)
{
    fmpz_t t, sq;
    fmpz_init(t);
    fmpz_init(sq);
    for (;;)
    {
        fmpz_randprime(n, state, bits, 0);
        if (modulus > 1 && fmpz_fdiv_ui(n, modulus) != residue)
            continue;
        fmpz_set_si(sq, D);
        fmpz_mod(sq, sq, n);
        if (fmpz_jacobi(sq, n) != 1)
            continue;
        fmpz_sqrtmod(sq, sq, n);
        if (!ecpp_cornacchia(t, v, n, D, sq))
            continue;
        if (vpar < 0 || (int) fmpz_is_even(v) == vpar)
            break;
    }
    fmpz_clear(t);
    fmpz_clear(sq);
}

/* checks that ecpp_class_poly_tower with the flags returns a root of H_Dt */
static void
_check_tower(slong D, slong Dt, int flags, const fmpz_t n, flint_rand_t state, const char * what)
{
    fmpz_mod_ctx_t ctx;
    fmpz_poly_t HZ;
    fmpz_mod_poly_t H;
    fmpz_t j, t;

    fmpz_init(j);
    fmpz_init(t);
    fmpz_mod_ctx_init(ctx, n);
    fmpz_poly_init(HZ);
    fmpz_mod_poly_init(H, ctx);

    if (!_ecpp_class_poly_tower(j, Dt, flags, state, ctx))
    {
        flint_printf("FAIL: no root found (%s), D = %wd, tower on %wd, flags %d\n", what, D, Dt, flags);
        flint_printf("n = "); fmpz_print(n); flint_printf("\n");
        flint_abort();
    }
    acb_modular_hilbert_class_poly(HZ, Dt);
    fmpz_mod_poly_set_fmpz_poly(H, HZ, ctx);
    fmpz_mod_poly_evaluate_fmpz(t, H, j, ctx);
    if (!fmpz_is_zero(t))
    {
        flint_printf("FAIL: not a root (%s), D = %wd, tower on %wd, flags %d\n", what, D, Dt, flags);
        flint_printf("n = "); fmpz_print(n); flint_printf("\n");
        flint_abort();
    }

    fmpz_mod_poly_clear(H, ctx);
    fmpz_poly_clear(HZ);
    fmpz_mod_ctx_clear(ctx);
    fmpz_clear(j);
    fmpz_clear(t);
}

TEST_FUNCTION_START(ecpp_class_poly_tower, state)
{
    slong iter, i;
    ulong primes[64];
    slong nprimes = 0;
    n_primes_t it;
    fmpz_t n, v;

    /* discriminants whose class number has a given prime factor p >= 5 at
       a level above the bottom, for the Kummer descents: h(-119) = 10,
       h(-215) = 14, h(-591) = 22 (odd D, hence j); h(-296) = 10,
       h(-404) = 14, h(-1256) = 26 (even D, Weber; -1256 also with j) */
    static const slong kummer_D[] = {-119, -215, -591, -296, -404, -1256, -1256};
    static const slong kummer_p[] = {5, 7, 11, 5, 7, 13, 13};
    /* even discriminants covering every class of m = -D/4 mod 8 in which
       Weber's function is used: m = 1, 2, 3, 5, 6, 7 */
    static const slong weber_D[] = {-4 * 33, -4 * 10, -4 * 35, -4 * 21, -4 * 22, -4 * 15,
                                    -4 * 105, -4 * 130, -4 * 195, -4 * 165, -4 * 462, -4 * 231};

    fmpz_init(n);
    fmpz_init(v);

    n_primes_init(it);
    n_primes_next(it);
    for (i = 0; i < 64; i++)
        primes[nprimes++] = n_primes_next(it);
    n_primes_clear(it);

    /* random discriminants from the pool, the default flags */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        ecpp_disc_struct * discs;
        slong ndiscs, k;
        const ecpp_disc_struct * d;

        ndiscs = ecpp_disc_table(&discs, primes, nprimes, 40000, 48, 48, 0.0);
        do
            k = n_randint(state, ndiscs);
        while (discs[k].h < 2);
        d = discs + k;

        _split_prime(n, v, d->D, 60 + n_randint(state, 100), 1, 0, -1, state);
        _check_tower(d->D, d->D, 0, n, state, "pool, default");

        /* odd D with v even: the order of conductor 2 with Weber's f */
        if (d->D % 2 != 0 && fmpz_is_even(v))
            _check_tower(d->D, 4 * d->D, 0, n, state, "conductor 2");

        flint_free(discs);
    }

    /* the Kummer descents: n = 1 mod p (in F_n) and n = -1 mod p (in
       F_{n^2}), and the powering when neither applies */
    for (i = 0; i < 7; i++)
    {
        slong D = kummer_D[i], p = kummer_p[i];
        int wflags = (D % 4 == 0 && i != 6) ? 0 : ECPP_TOWER_J;

        _split_prime(n, v, D, 100 + n_randint(state, 100), p, 1, -1, state);
        _check_tower(D, D, ECPP_TOWER_KUMMER | wflags, n, state, "Kummer in F_n");
        _split_prime(n, v, D, 100 + n_randint(state, 100), p, p - 1, -1, state);
        _check_tower(D, D, ECPP_TOWER_KUMMER | wflags, n, state, "Kummer in F_(n^2)");
        _split_prime(n, v, D, 100 + n_randint(state, 100), p, 2, -1, state);
        _check_tower(D, D, ECPP_TOWER_KUMMER | wflags, n, state, "Kummer data, powering");
    }

    /* class number one, and a Kummer descent needing the discrete
       logarithm in the p-Sylow subgroup (p^2 | n - 1) */
    _split_prime(n, v, -7, 60 + n_randint(state, 60), 1, 0, -1, state);
    _check_tower(-7, -7, 0, n, state, "h = 1");
    _split_prime(n, v, -296, 100 + n_randint(state, 60), 25, 1, -1, state);
    _check_tower(-296, -296, ECPP_TOWER_KUMMER, n, state, "Kummer, 25 | n - 1");
    _split_prime(n, v, -404, 100 + n_randint(state, 60), 49, 1, -1, state);
    _check_tower(-404, -404, ECPP_TOWER_KUMMER, n, state, "Kummer, 49 | n - 1");

    /* Weber's function in every residue class of m, against the j path */
    for (i = 0; i < 12; i++)
    {
        slong D = weber_D[i];
        _split_prime(n, v, D, 60 + n_randint(state, 100), 1, 0, -1, state);
        _check_tower(D, D, 0, n, state, "Weber");
        _check_tower(D, D, ECPP_TOWER_J, n, state, "j on an even D");
    }

    /* odd D of both residue classes mod 8 through the order of conductor
       2 (v even: always for D = 1 mod 8, half of the time for D = 5 mod
       8), and with j when v is odd (D = 5 mod 8 only) */
    for (iter = 0; iter < 4; iter++)
    {
        static const slong odd_D[] = {-1147, -455, -615, -1155};
        slong D = odd_D[iter];
        _split_prime(n, v, D, 60 + n_randint(state, 100), 1, 0, 1, state);
        _check_tower(D, 4 * D, 0, n, state, "conductor 2, v even");
        if ((-D) % 8 == 3)
        {
            _split_prime(n, v, D, 60 + n_randint(state, 100), 1, 0, 0, state);
            _check_tower(D, D, 0, n, state, "j, v odd");
        }
    }

    fmpz_clear(n);
    fmpz_clear(v);

    TEST_FUNCTION_END(state);
}
