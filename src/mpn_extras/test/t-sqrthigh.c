/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "mpn_extras.h"

/*
    flint_mpn_sqrthigh / rsqrthigh against the exact values computed with
    fmpz: the n + 1 limbs lie in [B^(n+1)/2, B^(n+1)) and within the error
    bound of S B^(n+1) resp. M B^(n+1); the _exact variants give the floor
    and the inexact flag.
*/

static slong
_sqrthigh_random_length(flint_rand_t state, slong big)
{
    switch (n_randint(state, 6))
    {
        case 0: return 1 + n_randint(state, big);
        case 1: return 1 + n_randint(state, 40);
        default: return 1 + n_randint(state, 8);
    }
}

static void
_sqrthigh_random_normalised(nn_ptr x, slong n, flint_rand_t state)
{
    slong i;

    if (n_randint(state, 4) == 0)
    {
        flint_mpn_zero(x, n);
        if (n_randint(state, 2))
            for (i = 0; i < n; i++)
                x[i] = n_randint(state, 3) ? UWORD_MAX : n_randtest(state);
    }
    else
        flint_mpn_rrandom(x, state, n);

    x[n - 1] |= (UWORD(1) << (FLINT_BITS - 1));
}

/* F = floor(S B^(n+1+extra)) = isqrt(floor(a B^(2(n+1+extra)) / (B^an 2^t)))
   (rsqrt = 0), or F = floor(M B^(n+1+extra)) =
   isqrt(floor(B^(2(n+1+extra)) B^an / (2^(2-t) a))) (rsqrt = 1) */
static void
exact_root(fmpz_t F, nn_srcptr a, slong an, int t, slong n, slong extra, int rsqrt)
{
    fmpz_t A, X;
    fmpz_init(A);
    fmpz_init(X);
    fmpz_set_ui_array(A, a, an);
    if (!rsqrt)
    {
        fmpz_mul_2exp(X, A, 2 * FLINT_BITS * (n + 1 + extra));
        fmpz_set_ui(A, 1);
        fmpz_mul_2exp(A, A, FLINT_BITS * an + t);
        fmpz_fdiv_q(X, X, A);
    }
    else
    {
        fmpz_mul_2exp(A, A, 2 - t);
        fmpz_set_ui(X, 1);
        fmpz_mul_2exp(X, X, FLINT_BITS * (2 * (n + 1 + extra) + an));
        fmpz_fdiv_q(X, X, A);
    }
    fmpz_sqrt(F, X);
    fmpz_clear(A);
    fmpz_clear(X);
}

/* whether S B^(n+1) resp. M B^(n+1) is an integer: the exact square equals
   the (exact rational) radicand */
static int
root_is_exact(const fmpz_t F, nn_srcptr a, slong an, int t, slong n, int rsqrt)
{
    fmpz_t A, X;
    int r;
    fmpz_init(A);
    fmpz_init(X);
    fmpz_set_ui_array(A, a, an);
    fmpz_mul(X, F, F);
    if (!rsqrt)
    {
        /* F^2 B^an 2^t == a B^(2(n+1)) */
        fmpz_mul_2exp(X, X, FLINT_BITS * an + t);
        fmpz_mul_2exp(A, A, 2 * FLINT_BITS * (n + 1));
    }
    else
    {
        /* F^2 2^(2-t) a == B^(2(n+1) + an) */
        fmpz_mul(X, X, A);
        fmpz_mul_2exp(X, X, 2 - t);
        fmpz_set_ui(A, 1);
        fmpz_mul_2exp(A, A, FLINT_BITS * (2 * (n + 1) + an));
    }
    r = fmpz_equal(X, A);
    fmpz_clear(A);
    fmpz_clear(X);
    return r;
}

TEST_FUNCTION_START(flint_mpn_sqrthigh, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        slong an, n, i, bound;
        nn_ptr a, q, qe;
        fmpz_t Q, F, D;
        int t, rsqrt, inexact;

        n = _sqrthigh_random_length(state, 300);
        an = _sqrthigh_random_length(state, 700);
        t = n_randint(state, 2);
        rsqrt = n_randint(state, 2);

        a = flint_malloc(sizeof(ulong) * an);
        q = flint_malloc(sizeof(ulong) * (n + 1));
        qe = flint_malloc(sizeof(ulong) * (n + 1));
        fmpz_init(Q);
        fmpz_init(F);
        fmpz_init(D);

        _sqrthigh_random_normalised(a, an, state);

        if (n_randint(state, 8) == 0)
        {
            /* a perfect square (of a fraction with a few limbs) */
            slong cn = 1 + n_randint(state, FLINT_MAX(1, FLINT_MIN(an / 2, n)));
            nn_ptr c = flint_malloc(sizeof(ulong) * cn), p = flint_malloc(sizeof(ulong) * 2 * cn);
            slong lead;
            flint_mpn_rrandom(c, state, cn);
            c[cn - 1] |= (UWORD(1) << (FLINT_BITS - 1));
            flint_mpn_sqr(p, c, cn);
            lead = flint_clz(p[2 * cn - 1]);
            if (lead)
                mpn_lshift(p, p, 2 * cn, lead);
            t = lead & 1;   /* the shift of the odd exponent */
            flint_mpn_zero(a, an);
            flint_mpn_copyi(a + an - FLINT_MIN(an, 2 * cn), p + 2 * cn - FLINT_MIN(an, 2 * cn), FLINT_MIN(an, 2 * cn));
            flint_free(c);
            flint_free(p);
        }

        else if (!rsqrt && n_randint(state, 8) == 0 && an == 2 * (n + 1) + 1)
        {
            /* a tie of the exact root: S B^(n+1) = F + 1/2 for a random
               F, i.e. a = (2 F + 1)^2 B / 2^(2 - t) with t = 1 unless
               (2 F + 1)^2 >= 2 B^(2n+2) */
            fmpz_t G;
            fmpz_init(G);
            fmpz_randbits(G, state, FLINT_BITS * (n + 1) - 1);
            fmpz_abs(G, G);
            fmpz_setbit(G, FLINT_BITS * (n + 1) - 1);
            fmpz_mul_2exp(G, G, 1);
            fmpz_add_ui(G, G, 1);
            fmpz_mul(G, G, G);
            t = (fmpz_bits(G) <= 2 * FLINT_BITS * (n + 1) + 1);
            fmpz_mul_2exp(G, G, FLINT_BITS - 2 + t);
            FLINT_ASSERT(fmpz_bits(G) == FLINT_BITS * an);
            fmpz_get_ui_array(a, an, G);
            fmpz_clear(G);
        }

        /* rsqrthigh requires M < 1 */
        if (rsqrt && t == 1 && a[an - 1] == (UWORD(1) << (FLINT_BITS - 1)) && flint_mpn_zero_p(a, an - 1))
            t = 0;

        for (i = 0; i <= n; i++)
            q[i] = qe[i] = UWORD(0xdeadbeef);

        if (rsqrt)
        {
            flint_mpn_rsqrthigh(q, a, an, t, n);
            bound = FLINT_MPN_RSQRTHIGH_GUARD_ERR_BOUND;
        }
        else
        {
            flint_mpn_sqrthigh(q, a, an, t, n);
            bound = FLINT_MPN_SQRTHIGH_GUARD_ERR_BOUND;
        }

        if (!(q[n] >> (FLINT_BITS - 1)))
            TEST_FUNCTION_FAIL("not normalised: an = %wd, n = %wd, t = %d, rsqrt = %d\n", an, n, t, rsqrt);

        exact_root(F, a, an, t, n, 1, rsqrt);
        fmpz_set_ui_array(Q, q, n + 1);
        fmpz_mul_2exp(Q, Q, FLINT_BITS);
        fmpz_sub(D, Q, F);
        fmpz_abs(D, D);
        fmpz_set_ui(Q, bound);
        fmpz_mul_2exp(Q, Q, FLINT_BITS);
        fmpz_add_ui(Q, Q, 1);
        if (fmpz_cmp(D, Q) > 0)
        {
            flint_printf("FAIL: error bound\n");
            flint_printf("an = %wd, n = %wd, t = %d, rsqrt = %d\n", an, n, t, rsqrt);
            flint_printf("|Q B - F| = "); fmpz_print(D); flint_printf("\n");
            flint_abort();
        }

        if (rsqrt)
            inexact = flint_mpn_rsqrthigh_exact(qe, a, an, t, n);
        else
            inexact = flint_mpn_sqrthigh_exact(qe, a, an, t, n);

        exact_root(F, a, an, t, n, 0, rsqrt);
        fmpz_set_ui_array(Q, qe, n + 1);
        if (!fmpz_equal(Q, F) || inexact != !root_is_exact(F, a, an, t, n, rsqrt))
        {
            flint_printf("FAIL: exact root\n");
            flint_printf("an = %wd, n = %wd, t = %d, rsqrt = %d, inexact = %d\n", an, n, t, rsqrt, inexact);
            flint_abort();
        }

        if (!rsqrt)
        {
            /* the code against the half point: (2 F + 1)^2 B^an 2^t
               compared with 4 a B^(2(n+1)) */
            int c, expected;

            for (i = 0; i <= n; i++)
                qe[i] = UWORD(0xdeadbeef);
            inexact = flint_mpn_sqrthigh_exact_half(qe, a, an, t, n);
            fmpz_set_ui_array(Q, qe, n + 1);
            if (!fmpz_equal(Q, F))
            {
                flint_printf("FAIL: exact root (half)\n");
                flint_printf("an = %wd, n = %wd, t = %d\n", an, n, t);
                flint_abort();
            }
            fmpz_mul_2exp(Q, F, 1);
            fmpz_add_ui(Q, Q, 1);
            fmpz_mul(Q, Q, Q);
            fmpz_mul_2exp(Q, Q, FLINT_BITS * an + t);
            fmpz_set_ui_array(D, a, an);
            fmpz_mul_2exp(D, D, 2 * FLINT_BITS * (n + 1) + 2);
            c = fmpz_cmp(Q, D);
            expected = (c < 0) ? 3 : (c == 0) ? 2 : (inexact != 0);
            if (inexact != expected)
            {
                flint_printf("FAIL: half code\n");
                flint_printf("an = %wd, n = %wd, t = %d, code = %d, expected = %d\n", an, n, t, inexact, expected);
                flint_abort();
            }
        }

        flint_free(a);
        flint_free(q);
        flint_free(qe);
        fmpz_clear(Q);
        fmpz_clear(F);
        fmpz_clear(D);
    }

    TEST_FUNCTION_END(state);
}
