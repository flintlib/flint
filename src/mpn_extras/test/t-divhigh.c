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
#include "mpn_extras/fp_small.h"

/*
    flint_mpn_divhigh / invhigh against the exact quotient computed with
    fmpz: the n + 1 limbs lie in [B^(n+1)/2, B^(n+1)) and within the
    error bound of 2^-s (A / Y) B^(n+1); flint_mpn_divhigh_exact gives the
    floor and the inexact flag.
*/

static slong
_divhigh_random_length(flint_rand_t state, slong big)
{
    switch (n_randint(state, 6))
    {
        case 0: return 1 + n_randint(state, big);
        case 1: return 1 + n_randint(state, 40);
        default: return 1 + n_randint(state, 8);
    }
}

static void
_divhigh_random_normalised(nn_ptr x, slong n, flint_rand_t state)
{
    slong i;

    if (n_randint(state, 4) == 0)
    {
        /* sparse: values near powers of two, equal top limbs */
        flint_mpn_zero(x, n);
        if (n_randint(state, 2))
            for (i = 0; i < n; i++)
                x[i] = n_randint(state, 3) ? UWORD_MAX : n_randtest(state);
    }
    else
        flint_mpn_rrandom(x, state, n);

    x[n - 1] |= (UWORD(1) << (FLINT_BITS - 1));
}

/* F = floor(2^-s (A / Y) B^(n+1+extra)) */
static void
exact_quotient(fmpz_t F, nn_srcptr a, slong an, nn_srcptr b, slong bn, ulong s, slong n, slong extra)
{
    fmpz_t A, Y;
    fmpz_init(A);
    fmpz_init(Y);
    fmpz_set_ui_array(A, a, an);
    fmpz_set_ui_array(Y, b, bn);
    /* A / Y = (a / B^an) / (b / B^bn): a B^(bn - an + n + 1 + extra) / b */
    fmpz_mul_2exp(A, A, FLINT_BITS * (bn + n + 1 + extra));
    fmpz_mul_2exp(Y, Y, FLINT_BITS * an + s);
    fmpz_fdiv_q(F, A, Y);
    fmpz_clear(A);
    fmpz_clear(Y);
}

TEST_FUNCTION_START(flint_mpn_divhigh, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        slong an, bn, n, i;
        nn_ptr a, b, q, qe;
        fmpz_t Q, F, D;
        ulong s, se;
        int inv, inexact, ok;

        n = _divhigh_random_length(state, 300);
        an = _divhigh_random_length(state, 320);
        bn = _divhigh_random_length(state, 320);
        inv = (n_randint(state, 6) == 0);

        a = flint_malloc(sizeof(ulong) * an);
        b = flint_malloc(sizeof(ulong) * bn);
        q = flint_malloc(sizeof(ulong) * (n + 1));
        qe = flint_malloc(sizeof(ulong) * (n + 1));
        fmpz_init(Q);
        fmpz_init(F);
        fmpz_init(D);

        _divhigh_random_normalised(a, an, state);
        _divhigh_random_normalised(b, bn, state);
        if (inv)
        {
            an = 1;
            a[0] = UWORD(1) << (FLINT_BITS - 1);
        }
        else if (n_randint(state, 8) == 0)
        {
            /* equal operands, or equal top limbs */
            slong m = FLINT_MIN(an, bn);
            flint_mpn_copyi(a + an - m, b + bn - m, m);
        }
        else if (n_randint(state, 8) == 0)
        {
            /* an exact quotient: a = b * c truncated to a fraction */
            slong cn = 1 + n_randint(state, 3);
            nn_ptr c = flint_malloc(sizeof(ulong) * cn), p;
            flint_mpn_rrandom(c, state, cn);
            c[cn - 1] |= 1;
            p = flint_malloc(sizeof(ulong) * (bn + cn));
            if (bn >= cn)
                flint_mpn_mul(p, b, bn, c, cn);
            else
                flint_mpn_mul(p, c, cn, b, bn);
            an = bn + cn - (p[bn + cn - 1] == 0);
            flint_free(a);
            a = flint_malloc(sizeof(ulong) * an);
            flint_mpn_copyi(a, p, an);
            {
                slong lead = flint_clz(a[an - 1]);
                if (lead)
                    mpn_lshift(a, a, an, lead);
            }
            flint_free(c);
            flint_free(p);
        }

        for (i = 0; i <= n; i++)
            q[i] = qe[i] = UWORD(0xdeadbeef);

        if (inv)
            s = flint_mpn_invhigh(q, b, bn, n);
        else
            s = flint_mpn_divhigh(q, a, an, b, bn, n);

        if (s != 0 && s != 1)
            TEST_FUNCTION_FAIL("shift %wu\n", s);

        if (!(q[n] >> (FLINT_BITS - 1)))
            TEST_FUNCTION_FAIL("not normalised: an = %wd, bn = %wd, n = %wd\n", an, bn, n);

        /* |Q B - F| <= bound B + 1 with F = floor(T B) */
        exact_quotient(F, a, an, b, bn, s, n, 1);
        fmpz_set_ui_array(Q, q, n + 1);
        fmpz_mul_2exp(Q, Q, FLINT_BITS);
        fmpz_sub(D, Q, F);
        fmpz_abs(D, D);
        fmpz_set_ui(Q, FLINT_MPN_DIVHIGH_GUARD_ERR_BOUND);
        fmpz_mul_2exp(Q, Q, FLINT_BITS);
        fmpz_add_ui(Q, Q, 1);
        if (fmpz_cmp(D, Q) > 0)
        {
            flint_printf("FAIL: error bound\n");
            flint_printf("an = %wd, bn = %wd, n = %wd, s = %wu, inv = %d\n", an, bn, n, s, inv);
            flint_printf("|Q B - F| = "); fmpz_print(D); flint_printf("\n");
            flint_abort();
        }

        /* the exact floor */
        inexact = flint_mpn_divhigh_exact(qe, &se, a, an, b, bn, n);
        exact_quotient(F, a, an, b, bn, se, n, 0);
        fmpz_set_ui_array(Q, qe, n + 1);
        ok = fmpz_equal(Q, F);
        if (ok)
        {
            /* exact iff 2^s b B^an divides a B^(bn + n + 1) */
            fmpz_t A, Y;
            fmpz_init(A);
            fmpz_init(Y);
            fmpz_set_ui_array(A, a, an);
            fmpz_mul_2exp(A, A, FLINT_BITS * (bn + n + 1));
            fmpz_set_ui_array(Y, b, bn);
            fmpz_mul_2exp(Y, Y, FLINT_BITS * an + se);
            fmpz_fdiv_r(D, A, Y);
            ok = (inexact == !fmpz_is_zero(D));
            fmpz_clear(A);
            fmpz_clear(Y);
        }
        if (!ok)
        {
            flint_printf("FAIL: exact quotient\n");
            flint_printf("an = %wd, bn = %wd, n = %wd, s = %wu, inexact = %d\n", an, bn, n, se, inexact);
            flint_printf("Q = "); fmpz_set_ui_array(Q, qe, n + 1); fmpz_print(Q); flint_printf("\n");
            flint_printf("F = "); fmpz_print(F); flint_printf("\n");
            flint_printf("a = "); for (i = an - 1; i >= 0; i--) flint_printf("%wx ", a[i]); flint_printf("\n");
            flint_printf("b = "); for (i = bn - 1; i >= 0; i--) flint_printf("%wx ", b[i]); flint_printf("\n");
            flint_abort();
        }

        flint_free(a);
        flint_free(b);
        flint_free(q);
        flint_free(qe);
        fmpz_clear(Q);
        fmpz_clear(F);
        fmpz_clear(D);
    }

    /* the register kernel for a three-limb divisor against the exact
       floor (flint_mpn_divhigh_exact, with n = 2) */
    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        ulong a[2], b[3], q[3], qe[3], q2, q1, q0, s, se;
        slong an;
        int inexact, inexact_e;

        an = 1 + n_randint(state, 2);
        _divhigh_random_normalised(a + 2 - an, an, state);
        if (an == 1)
            a[0] = 0;
        _divhigh_random_normalised(b, 3, state);
        if (n_randint(state, 4) == 0)
        {
            /* dividends near the divisor */
            a[1] = b[2];
            if (n_randint(state, 2))
                a[0] = b[1];
        }

        inexact = _flint_mpn_fdiv_3_kernel(&q2, &q1, &q0, &s, a[1], a[0], b[2], b[1], b[0]);
        q[2] = q2; q[1] = q1; q[0] = q0;
        inexact_e = flint_mpn_divhigh_exact(qe, &se, a, 2, b, 3, 2);

        if (s != se || inexact != inexact_e || q[0] != qe[0] || q[1] != qe[1] || q[2] != qe[2])
        {
            flint_printf("FAIL: _flint_mpn_fdiv_3_kernel\n");
            flint_printf("a = %wx %wx\nb = %wx %wx %wx\n", a[1], a[0], b[2], b[1], b[0]);
            flint_printf("s = %wu (%wu), inexact = %d (%d)\n", s, se, inexact, inexact_e);
            flint_printf("q  = %wx %wx %wx\nqe = %wx %wx %wx\n", q[2], q[1], q[0], qe[2], qe[1], qe[0]);
            flint_abort();
        }
    }

    TEST_FUNCTION_END(state);
}
