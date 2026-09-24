/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "t-helpers.h"

/*
    y^2 = x^3 + a x + b over Z/n through a point (x, y) chosen first, so
    that b makes the point lie on the curve whatever n is.
*/
static int
build(gr_ctx_t R, gr_ec_ctx_t E, gr_ec_jac_point_t P, const fmpz_t n,
        const fmpz_t a, const fmpz_t x, const fmpz_t y)
{
    gr_ptr g4, g6, gx, gy;
    int ok = 0;

    gr_ctx_init_fmpz_mod(R, n);

    GR_TMP_INIT4(g4, g6, gx, gy, R);

    if (gr_set_fmpz(g4, a, R) == GR_SUCCESS
            && gr_set_fmpz(gx, x, R) == GR_SUCCESS
            && gr_set_fmpz(gy, y, R) == GR_SUCCESS)
    {
        /* b = y^2 - x^3 - a x */
        int st = GR_SUCCESS;
        gr_ptr t = g6;

        st |= gr_sqr(t, gy, R);
        st |= gr_sqr(gy, gx, R);
        st |= gr_mul(gy, gy, gx, R);
        st |= gr_sub(t, t, gy, R);
        st |= gr_mul(gy, g4, gx, R);
        st |= gr_sub(g6, t, gy, R);

        if (st == GR_SUCCESS
                && gr_ec_ctx_init_short_weierstrass(E, R, g4, g6) == GR_SUCCESS)
            ok = 1;
    }

    GR_TMP_CLEAR4(g4, g6, gx, gy, R);

    if (ok)
    {
        gr_ec_jac_point_init(P, E);

        if (gr_set_fmpz(GR_EC_JAC_POINT_X(P, E), x, R) != GR_SUCCESS
                || gr_set_fmpz(GR_EC_JAC_POINT_Y(P, E), y, R) != GR_SUCCESS
                || gr_one(GR_EC_JAC_POINT_Z(P, E), R) != GR_SUCCESS)
            ok = 0;

        P->is_infinity = T_FALSE;
    }

    return ok;
}

/*
    Over a prime modulus every nonzero element is a unit, so every branch
    the group law takes really is entitled to its decision and the witness
    must always come out a unit. Anything else would be a false alarm, and
    a false alarm makes a primality prover give up on a prime.
*/
static void
check_sound_over_a_prime(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        fmpz_t p, a, x, y, k, g, wf;
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ec_jac_point_t P, Q;
        gr_ptr w;

        fmpz_init(p); fmpz_init(a); fmpz_init(x); fmpz_init(y);
        fmpz_init(k); fmpz_init(g); fmpz_init(wf);

        fmpz_randprime(p, state, 10 + n_randint(state, 54), 0);
        fmpz_randm(a, state, p);
        fmpz_randm(x, state, p);
        fmpz_randm(y, state, p);
        fmpz_randm(k, state, p);

        if (build(R, E, P, p, a, x, y))
        {
            gr_ec_jac_point_init(Q, E);
            w = gr_heap_init(R);

            FLINT_TEST(gr_one(w, R) == GR_SUCCESS);

            if (gr_ec_jac_point_mul_fmpz_witness(Q, w, P, k, E) == GR_SUCCESS)
            {
                gr_ec_jac_point_t Q2;

                FLINT_TEST(gr_get_fmpz(wf, w, R) == GR_SUCCESS);
                fmpz_gcd(g, wf, p);
                FLINT_TEST(fmpz_is_one(g));

                /* and the point is the same one the plain ladder gives */
                gr_ec_jac_point_init(Q2, E);

                if (gr_ec_jac_point_mul_fmpz(Q2, P, k, E) == GR_SUCCESS)
                    FLINT_TEST(gr_ec_jac_point_equal(Q, Q2, E) != T_FALSE);

                gr_ec_jac_point_clear(Q2, E);
            }

            gr_heap_clear(w, R);
            gr_ec_jac_point_clear(Q, E);
            gr_ec_ctx_clear(E);
        }

        gr_ctx_clear(R);
        fmpz_clear(p); fmpz_clear(a); fmpz_clear(x); fmpz_clear(y);
        fmpz_clear(k); fmpz_clear(g); fmpz_clear(wf);
    }
}

/*
    Over a composite modulus the witness is supposed to catch the group law
    out, and when it does the gcd is a factor. This is elliptic curve
    factorization: multiply by the prime powers up to a bound and wait for
    a branch to be taken on a quantity that is zero modulo one factor of n
    and not the other.

    Not every curve finds a factor within the bound, which is the nature of
    the method; what is checked is that whenever the witness does report
    something, it reports a genuine proper factor.
*/
static void
check_detects_over_a_composite(flint_rand_t state)
{
    slong iter, found = 0;

    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        fmpz_t p, q, n, a, x, y, k, g, wf;
        gr_ctx_t R;
        gr_ec_ctx_t E;
        gr_ec_jac_point_t P, T;
        gr_ptr w;
        ulong pr;
        int stop = 0;

        fmpz_init(p); fmpz_init(q); fmpz_init(n); fmpz_init(a);
        fmpz_init(x); fmpz_init(y); fmpz_init(k); fmpz_init(g); fmpz_init(wf);

        fmpz_set_ui(p, n_randprime(state, 16 + n_randint(state, 6), 1));
        fmpz_set_ui(q, n_randprime(state, 16 + n_randint(state, 6), 1));
        fmpz_mul(n, p, q);

        if (fmpz_equal(p, q))
            goto next;

        fmpz_randm(a, state, n);
        fmpz_randm(x, state, n);
        fmpz_randm(y, state, n);

        if (!build(R, E, P, n, a, x, y))
            goto next;

        gr_ec_jac_point_init(T, E);
        w = gr_heap_init(R);
        FLINT_TEST(gr_one(w, R) == GR_SUCCESS);

        for (pr = 2; pr < 2000 && !stop; pr = n_nextprime(pr, 1))
        {
            ulong e = pr;

            while (e < 2000)
                e *= pr;

            fmpz_set_ui(k, e / pr);

            if (gr_ec_jac_point_mul_fmpz_witness(T, w, P, k, E) != GR_SUCCESS)
                stop = 1;
            else
                FLINT_TEST(gr_ec_jac_point_set(P, T, E) == GR_SUCCESS);

            FLINT_TEST(gr_get_fmpz(wf, w, R) == GR_SUCCESS);
            fmpz_gcd(g, wf, n);

            if (!fmpz_is_one(g))
            {
                /* whatever it says has to be true */
                FLINT_TEST(fmpz_divisible(n, g));

                if (!fmpz_equal(g, n))
                {
                    FLINT_TEST(fmpz_equal(g, p) || fmpz_equal(g, q));
                    found++;
                }

                stop = 1;
            }
        }

        gr_heap_clear(w, R);
        gr_ec_jac_point_clear(T, E);
        gr_ec_ctx_clear(E);
        gr_ctx_clear(R);

next:
        fmpz_clear(p); fmpz_clear(q); fmpz_clear(n); fmpz_clear(a);
        fmpz_clear(x); fmpz_clear(y); fmpz_clear(k); fmpz_clear(g); fmpz_clear(wf);
    }

    /* the method is not certain, but over this many curves it does work */
    FLINT_TEST(found > 0);
}

TEST_FUNCTION_START(gr_ec_witness, state)
{
    check_sound_over_a_prime(state);
    check_detects_over_a_composite(state);

    TEST_FUNCTION_END(state);
}
