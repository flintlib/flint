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
#include "t-helpers.h"

/* |N - (q + 1)| <= 2 sqrt(q), which is Hasse's theorem */
static void
check_hasse(const fmpz_t N, const fmpz_t q)
{
    fmpz_t t, bound;

    fmpz_init(t);
    fmpz_init(bound);

    fmpz_sub(t, N, q);
    fmpz_sub_ui(t, t, 1);           /* -trace */
    fmpz_mul(t, t, t);
    fmpz_mul_ui(bound, q, 4);

    FLINT_TEST(fmpz_cmp(t, bound) <= 0);

    fmpz_clear(t);
    fmpz_clear(bound);
}

/* the cardinality of the base ring, through the new generic entry point */
static void
check_base_ring_cardinality(void)
{
    gr_ctx_t R;
    fmpz_t c;

    fmpz_init(c);

    /* Z/n, prime or not, and without needing a primality test */
    GR_MUST_SUCCEED(gr_ctx_init_nmod(R, 97));
    FLINT_TEST(gr_ctx_cardinality_fmpz(c, R) == GR_SUCCESS);
    FLINT_TEST(fmpz_equal_ui(c, 97));
    gr_ctx_clear(R);

    GR_MUST_SUCCEED(gr_ctx_init_nmod(R, 100));
    FLINT_TEST(gr_ctx_cardinality_fmpz(c, R) == GR_SUCCESS);
    FLINT_TEST(fmpz_equal_ui(c, 100));
    gr_ctx_clear(R);

    /* a finite field of prime power order */
    gr_ctx_init_fq_nmod(R, 7, 3, "a");
    FLINT_TEST(gr_ctx_cardinality_fmpz(c, R) == GR_SUCCESS);
    FLINT_TEST(fmpz_equal_ui(c, 343));
    gr_ctx_clear(R);

    gr_ctx_init_fq_nmod(R, 2, 5, "a");
    FLINT_TEST(gr_ctx_cardinality_fmpz(c, R) == GR_SUCCESS);
    FLINT_TEST(fmpz_equal_ui(c, 32));
    gr_ctx_clear(R);

    /* a large modulus: no primality test is needed to know the size */
    {
        fmpz_t p;
        fmpz_init(p);
        fmpz_set_str(p, "115792089237316195423570985008687907853269984665640564039457584007913129640233", 10);
        gr_ctx_init_fmpz_mod(R, p);
        FLINT_TEST(gr_ctx_cardinality_fmpz(c, R) == GR_SUCCESS);
        FLINT_TEST(fmpz_equal(c, p));
        gr_ctx_clear(R);

        FLINT_TEST(gr_ctx_init_mpn_mod(R, p) == GR_SUCCESS);
        FLINT_TEST(gr_ctx_cardinality_fmpz(c, R) == GR_SUCCESS);
        FLINT_TEST(fmpz_equal(c, p));
        gr_ctx_clear(R);
        fmpz_clear(p);
    }

    /* an infinite structure has no cardinality, which is a domain answer */
    gr_ctx_init_fmpz(R);
    FLINT_TEST(gr_ctx_cardinality_fmpz(c, R) == GR_DOMAIN);
    gr_ctx_clear(R);

    gr_ctx_init_fmpq(R);
    FLINT_TEST(gr_ctx_cardinality_fmpz(c, R) == GR_DOMAIN);
    gr_ctx_clear(R);

    fmpz_clear(c);
}

/*
    The three algorithms must agree with each other, satisfy Hasse, kill
    every point of the curve, and be what the generic interface reports.
*/
static void
check_algorithms_agree(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 120 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        fmpz_t q, nn, nb, ns, ng;
        slong t;
        ulong pp;
        int deg = 1;

        /* small fields, including characteristic 2 and 3 */
        switch (n_randint(state, 4))
        {
            case 0: pp = n_randprime(state, 5 + n_randint(state, 7), 1); break;
            case 1: pp = n_randprime(state, 4 + n_randint(state, 4), 1); break;
            case 2: pp = 2; deg = 2 + n_randint(state, 5); break;
            default: pp = 3; deg = 1 + n_randint(state, 4); break;
        }

        if (deg == 1)
        {
            if (gr_ctx_init_nmod(R, pp) != GR_SUCCESS)
                continue;
        }
        else
            gr_ctx_init_fq_nmod(R, pp, deg, "a");

        fmpz_init(q); fmpz_init(nn); fmpz_init(nb); fmpz_init(ns); fmpz_init(ng);

        if (gr_ec_ctx_init_randtest(E, state, R) != GR_SUCCESS)
            goto next;

        FLINT_TEST(gr_ctx_cardinality_fmpz(q, R) == GR_SUCCESS);

        /* naive is the reference: it never touches the group law */
        if (gr_ec_ctx_cardinality_naive(nn, E) != GR_SUCCESS)
            goto next_curve;

        check_hasse(nn, q);

        /* BSGS may legitimately fail to disambiguate over a tiny field */
        if (gr_ec_ctx_cardinality_bsgs(nb, E) == GR_SUCCESS)
            FLINT_TEST(fmpz_equal(nn, nb));

        /* Schoof only claims the short model away from 2 and 3 */
        if (gr_ec_ctx_cardinality_schoof(ns, E) == GR_SUCCESS)
            FLINT_TEST(fmpz_equal(nn, ns));

        /* the CM path declines unless the j-invariant proves it applies */
        {
            fmpz_t nc;
            fmpz_init(nc);
            if (gr_ec_ctx_cardinality_cm(nc, E) == GR_SUCCESS)
                FLINT_TEST(fmpz_equal(nn, nc));
            fmpz_clear(nc);
        }

        /* and the same number must come out of the generic interface */
        FLINT_TEST(gr_ctx_cardinality_fmpz(ng, E) == GR_SUCCESS);
        FLINT_TEST(fmpz_equal(nn, ng));

        /* N annihilates the group */
        for (t = 0; t < 6; t++)
        {
            gr_ec_point_t P, Q;
            gr_ec_point_init(P, E);
            gr_ec_point_init(Q, E);

            if (gr_ec_point_randtest(P, state, E) == GR_SUCCESS
                    && gr_ec_point_mul_fmpz(Q, P, nn, E) == GR_SUCCESS)
                FLINT_TEST(gr_ec_point_is_inf(Q, E) == T_TRUE);

            gr_ec_point_clear(Q, E);
            gr_ec_point_clear(P, E);
        }

next_curve:
        gr_ec_ctx_clear(E);
next:
        fmpz_clear(q); fmpz_clear(nn); fmpz_clear(nb); fmpz_clear(ns); fmpz_clear(ng);
        gr_ctx_clear(R);
    }
}

/*
    Two families whose answer is known in advance. Over F_p with p > 3,
    y^2 = x^3 + b is supersingular when p = 2 mod 3, and y^2 = x^3 + a x is
    supersingular when p = 3 mod 4; either way #E = p + 1 exactly.
*/
static void
check_supersingular(flint_rand_t state)
{
    slong iter;

    for (iter = 0; iter < 40 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t R;
        gr_ec_ctx_t E;
        fmpz_t n;
        gr_ptr a, z;
        ulong pp;
        int j0 = n_randint(state, 2);
        slong tries;

        /* p = 2 mod 3 for the j = 0 family, p = 3 mod 4 for j = 1728 */
        for (tries = 0; tries < 200; tries++)
        {
            pp = n_randprime(state, 6 + n_randint(state, 12), 1);
            if (pp > 3 && ((j0 && pp % 3 == 2) || (!j0 && pp % 4 == 3)))
                break;
        }

        if (pp <= 3 || gr_ctx_init_nmod(R, pp) != GR_SUCCESS)
            continue;

        fmpz_init(n);
        GR_TMP_INIT2(a, z, R);

        if (gr_randtest_not_zero(a, state, R) != GR_SUCCESS
                || gr_zero(z, R) != GR_SUCCESS)
            goto next;

        /* y^2 = x^3 + a  or  y^2 = x^3 + a x */
        if (gr_ec_ctx_init_short_weierstrass(E, R, j0 ? z : a, j0 ? a : z)
                != GR_SUCCESS)
            goto next;

        if (gr_ec_ctx_cardinality_naive(n, E) == GR_SUCCESS)
            FLINT_TEST(fmpz_equal_ui(n, pp + 1));

        if (gr_ec_ctx_cardinality_bsgs(n, E) == GR_SUCCESS)
            FLINT_TEST(fmpz_equal_ui(n, pp + 1));

        if (gr_ec_ctx_cardinality_schoof(n, E) == GR_SUCCESS)
            FLINT_TEST(fmpz_equal_ui(n, pp + 1));

        /* these are exactly the curves the CM path is meant to recognise */
        FLINT_TEST(gr_ec_ctx_cardinality_cm(n, E) == GR_SUCCESS);
        FLINT_TEST(fmpz_equal_ui(n, pp + 1));

        gr_ec_ctx_clear(E);
next:
        GR_TMP_CLEAR2(a, z, R);
        fmpz_clear(n);
        gr_ctx_clear(R);
    }
}

/*
    Every class number one discriminant, against baby-step giant-step.

    The CM path reads the trace off a Cornacchia solution together with a
    handful of Jacobi symbols whose constants are tabulated per
    discriminant, so each entry of that table wants exercising: primes that
    split (ordinary, trace of either sign, reached through a random twist)
    and primes that are inert (supersingular).
*/
static const slong cm_test_j[] = {
    WORD(0), WORD(1728), WORD(-3375), WORD(8000), WORD(-32768),
    WORD(54000), WORD(287496), WORD(-884736), WORD(-12288000),
    WORD(16581375), WORD(-884736000), WORD(-147197952000),
    WORD(-262537412640768000)
};

static void
check_cm_discriminants(flint_rand_t state)
{
    slong idx, iter;

    for (idx = 0; idx < (slong) (sizeof(cm_test_j) / sizeof(slong)); idx++)
    {
        for (iter = 0; iter < 8 * flint_test_multiplier(); iter++)
        {
            gr_ctx_t R;
            gr_ec_ctx_t E;
            gr_ptr g4, g6;
            fmpz_t nc, nb;
            ulong pp, pinv, jm, a4, a6, c, c2, k;
            int built = 0;

            pp = n_randprime(state, 18 + n_randint(state, 7), 1);

            if (pp <= 2000 || gr_ctx_init_nmod(R, pp) != GR_SUCCESS)
                continue;

            pinv = n_preinvert_limb(pp);

            {
                fmpz_t fj, fp;
                fmpz_init_set_si(fj, cm_test_j[idx]);
                fmpz_init_set_ui(fp, pp);
                fmpz_mod(fj, fj, fp);
                jm = fmpz_get_ui(fj);
                fmpz_clear(fj);
                fmpz_clear(fp);
            }

            if (jm == 0)                            /* sextic twists */
            {
                a4 = 0;
                a6 = 1 + n_randint(state, pp - 1);
            }
            else if (jm == 1728 % pp)               /* quartic twists */
            {
                a4 = 1 + n_randint(state, pp - 1);
                a6 = 0;
            }
            else
            {
                /* a4 = 3 j (1728 - j), a6 = 2 j (1728 - j)^2, then a
                   random quadratic twist */
                k = n_submod(1728 % pp, jm, pp);
                a4 = n_mulmod2_preinv(n_mulmod2_preinv(3, jm, pp, pinv), k, pp, pinv);
                a6 = n_mulmod2_preinv(n_mulmod2_preinv(2, jm, pp, pinv),
                        n_mulmod2_preinv(k, k, pp, pinv), pp, pinv);

                c = 1 + n_randint(state, pp - 1);
                c2 = n_mulmod2_preinv(c, c, pp, pinv);
                a4 = n_mulmod2_preinv(a4, c2, pp, pinv);
                a6 = n_mulmod2_preinv(a6, n_mulmod2_preinv(c2, c, pp, pinv), pp, pinv);

                if (a4 == 0 || a6 == 0)
                {
                    gr_ctx_clear(R);
                    continue;
                }
            }

            GR_TMP_INIT2(g4, g6, R);

            if (gr_set_ui(g4, a4, R) == GR_SUCCESS
                    && gr_set_ui(g6, a6, R) == GR_SUCCESS
                    && gr_ec_ctx_init_short_weierstrass(E, R, g4, g6) == GR_SUCCESS)
                built = 1;

            GR_TMP_CLEAR2(g4, g6, R);

            if (!built)
            {
                gr_ctx_clear(R);
                continue;
            }

            fmpz_init(nc);
            fmpz_init(nb);

            /* the j-invariant proves the complex multiplication, so this
               one is not allowed to decline */
            FLINT_TEST(gr_ec_ctx_cardinality_cm(nc, E) == GR_SUCCESS);

            if (gr_ec_ctx_cardinality_bsgs(nb, E) == GR_SUCCESS)
                FLINT_TEST(fmpz_equal(nc, nb));

            fmpz_clear(nc);
            fmpz_clear(nb);
            gr_ec_ctx_clear(E);
            gr_ctx_clear(R);
        }
    }
}

/*
    Counting over F_{p^n} for n = 2 to 5.

    Nothing in baby-step giant-step or in Schoof is tied to a prime field,
    and the dispatcher should reach both of them over an extension; what is
    tied to a prime field is the complex multiplication shortcut, whose
    class number one theory is written for j in F_p, and it has to decline
    rather than answer. The naive walk is only affordable for the smallest
    of these.
*/
static void
check_extension_fields(flint_rand_t state)
{
    slong deg, iter;

    for (deg = 2; deg <= 5; deg++)
    {
        for (iter = 0; iter < 3 * flint_test_multiplier(); iter++)
        {
            gr_ctx_t R;
            gr_ec_ctx_t E;
            gr_ptr a4, a6;
            fmpz_t q, nn, nb, ns, nc, ng;
            slong tries;
            ulong p;
            int built = 0, have_n = 0, have_b = 0, have_s = 0;

            /* q = p^deg has to stay small: baby-step giant-step costs
               q^(1/4) group operations and Schoof works in a quotient of
               degree growing with log q, so a big p here would turn a unit
               test into a benchmark */
            p = n_nextprime(4 + n_randint(state, 40), 1);

            gr_ctx_init_fq_nmod(R, p, deg, "a");
            GR_TMP_INIT2(a4, a6, R);

            for (tries = 0; tries < 20 && !built; tries++)
                if (gr_randtest(a4, state, R) == GR_SUCCESS
                        && gr_randtest(a6, state, R) == GR_SUCCESS
                        && gr_ec_ctx_init_short_weierstrass(E, R, a4, a6) == GR_SUCCESS)
                    built = 1;

            GR_TMP_CLEAR2(a4, a6, R);

            if (!built)
            {
                gr_ctx_clear(R);
                continue;
            }

            fmpz_init(q); fmpz_init(nn); fmpz_init(nb);
            fmpz_init(ns); fmpz_init(nc); fmpz_init(ng);

            FLINT_TEST(gr_ctx_cardinality_fmpz(q, R) == GR_SUCCESS);

            /* the walk is O(q), and q is in the millions here, so only
               ask for it where it is genuinely cheap */
            have_n = (fmpz_cmp_ui(q, 100000) <= 0)
                && (gr_ec_ctx_cardinality_naive(nn, E) == GR_SUCCESS);
            have_b = (gr_ec_ctx_cardinality_bsgs(nb, E) == GR_SUCCESS);
            have_s = (gr_ec_ctx_cardinality_schoof(ns, E) == GR_SUCCESS);

            /* over an extension at least one of the two general
               algorithms has to deliver */
            FLINT_TEST(have_b || have_s);

            if (have_n) check_hasse(nn, q);
            if (have_b) check_hasse(nb, q);
            if (have_s) check_hasse(ns, q);

            if (have_n && have_b) FLINT_TEST(fmpz_equal(nn, nb));
            if (have_n && have_s) FLINT_TEST(fmpz_equal(nn, ns));
            if (have_b && have_s) FLINT_TEST(fmpz_equal(nb, ns));

            /* the CM path is for prime fields and must say so */
            FLINT_TEST(gr_ec_ctx_cardinality_cm(nc, E) == GR_UNABLE);

            /* and the generic entry point agrees with whichever ran */
            FLINT_TEST(gr_ctx_cardinality_fmpz(ng, E) == GR_SUCCESS);
            if (have_b) FLINT_TEST(fmpz_equal(ng, nb));
            if (have_s) FLINT_TEST(fmpz_equal(ng, ns));

            fmpz_clear(q); fmpz_clear(nn); fmpz_clear(nb);
            fmpz_clear(ns); fmpz_clear(nc); fmpz_clear(ng);
            gr_ec_ctx_clear(E);
            gr_ctx_clear(R);
        }
    }
}

/*
    A curve over F_{p^n} whose coefficients lie in the prime field is a
    base change, and the descent has to get the same answer as counting it
    where it stands. A curve that is not one has to be declined.
*/
static void
check_subfield(flint_rand_t state)
{
    slong deg, iter;

    for (deg = 2; deg <= 5; deg++)
    {
        for (iter = 0; iter < 4 * flint_test_multiplier(); iter++)
        {
            gr_ctx_t R;
            gr_ec_ctx_t E;
            gr_ptr a4, a6;
            fmpz_t ns, nb, dummy;
            slong tries;
            ulong p = n_nextprime(4 + n_randint(state, 60), 1);
            int built = 0;

            gr_ctx_init_fq_nmod(R, p, deg, "a");
            GR_TMP_INIT2(a4, a6, R);

            /* coefficients drawn from the prime field */
            for (tries = 0; tries < 30 && !built; tries++)
                if (gr_set_ui(a4, n_randint(state, p), R) == GR_SUCCESS
                        && gr_set_ui(a6, n_randint(state, p), R) == GR_SUCCESS
                        && gr_ec_ctx_init_short_weierstrass(E, R, a4, a6) == GR_SUCCESS)
                    built = 1;

            GR_TMP_CLEAR2(a4, a6, R);

            if (!built)
            {
                gr_ctx_clear(R);
                continue;
            }

            fmpz_init(ns);
            fmpz_init(nb);

            FLINT_TEST(gr_ec_ctx_cardinality_subfield(ns, E) == GR_SUCCESS);

            if (gr_ec_ctx_cardinality_bsgs(nb, E) == GR_SUCCESS)
                FLINT_TEST(fmpz_equal(ns, nb));

            /* and it must be what the dispatcher reports */
            FLINT_TEST(gr_ctx_cardinality_fmpz(nb, E) == GR_SUCCESS);
            FLINT_TEST(fmpz_equal(ns, nb));

            fmpz_clear(ns);
            fmpz_clear(nb);
            gr_ec_ctx_clear(E);

            /* over a prime field there is nothing to descend to */
            fmpz_init(dummy);
            gr_ctx_clear(R);

            if (gr_ctx_init_nmod(R, p) == GR_SUCCESS)
            {
                if (gr_ec_ctx_init_randtest(E, state, R) == GR_SUCCESS)
                {
                    FLINT_TEST(gr_ec_ctx_cardinality_subfield(dummy, E) == GR_UNABLE);
                    gr_ec_ctx_clear(E);
                }

                gr_ctx_clear(R);
            }

            fmpz_clear(dummy);
        }
    }
}

TEST_FUNCTION_START(gr_ec_cardinality, state)
{
    check_base_ring_cardinality();
    check_algorithms_agree(state);
    check_supersingular(state);
    check_cm_discriminants(state);
    check_extension_fields(state);
    check_subfield(state);

    TEST_FUNCTION_END(state);
}
