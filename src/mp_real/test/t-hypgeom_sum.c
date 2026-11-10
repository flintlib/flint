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
#include "fmpq.h"
#include "arb.h"
#include "mp_real.h"

/* mp_real_hypgeom_sum against the exact rational sum: random P, R and Q
   with positive coefficients (so that Q(k) > 0), convergent or not
   (deg R > deg Q), sometimes Q a large constant (the content-power mode
   with exact = 0); T/Q must contain the sum of the N' terms reported,
   N' = N for exact = 1 and N <= N' < N (1 + 1/12) otherwise, with an
   error below about 2^-(FLINT_BITS (n - 1)) times the largest term (as
   for mp_real_hypgeom_series, about n - 1 limbs). */

static void
_hs_set_int(mp_real_hypgeom_int_struct * x, const fmpz_t c)
{
    if (!COEFF_IS_MPZ(*c))
    {
        ulong * d = flint_malloc(sizeof(ulong));
        d[0] = FLINT_UABS(*c);
        x->d = d;
        x->n = (d[0] != 0);
        x->neg = (*c < 0);
    }
    else
    {
        mpz_ptr z = COEFF_TO_PTR(*c);
        ulong * d = flint_malloc(sizeof(ulong) * FLINT_ABS(z->_mp_size));
        flint_mpn_copyi(d, z->_mp_d, FLINT_ABS(z->_mp_size));
        x->d = d;
        x->n = FLINT_ABS(z->_mp_size);
        x->neg = (z->_mp_size < 0);
    }
}

static void
_hs_eval(fmpz_t v, const fmpz * f, slong len, slong k)
{
    slong i;
    fmpz_zero(v);
    for (i = len - 1; i >= 0; i--)
    {
        fmpz_mul_si(v, v, k);
        fmpz_add(v, v, f + i);
    }
}

TEST_FUNCTION_START(mp_real_hypgeom_sum, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        slong Plen, Qlen, Rlen, N, N2, n, k, i;
        int exact;
        fmpz P[4], Q[4], R[4];
        mp_real_hypgeom_int_struct Pi[4], Qi[4], Ri[4];
        fmpz_t v, w, U, V, Ar;
        fmpq_t S;
        slong Mexp;
        mp_real_t T, D;
        arb_t a, b;

        Plen = 1 + n_randint(state, 3);
        Qlen = 1 + n_randint(state, 3);
        Rlen = n_randint(state, 4);
        N = n_randint(state, 200);
        n = 1 + n_randint(state, 30);
        exact = n_randint(state, 2);

        for (i = 0; i < 4; i++)
        {
            fmpz_init(P + i);
            fmpz_init(Q + i);
            fmpz_init(R + i);
        }
        for (i = 0; i < Plen; i++)
            fmpz_randtest(P + i, state, 1 + n_randint(state, 100));
        for (i = 0; i < Qlen; i++)
        {
            fmpz_randtest_unsigned(Q + i, state, 1 + n_randint(state, 100));
            fmpz_add_ui(Q + i, Q + i, 1);
        }
        for (i = 0; i < Rlen; i++)
            fmpz_randtest(R + i, state, 1 + n_randint(state, 60));
        if (iter % 8 == 0)
        {
            /* a large constant Q, small R, many terms at a high
               precision: the content-power mode */
            N = 100 + n_randint(state, 1200);
            n = 64 + n_randint(state, 60);
            Plen = 1;
            fmpz_randtest(P + 0, state, 10);
            Qlen = 1;
            fmpz_randbits(Q + 0, state, 200 + n_randint(state, 200));
            fmpz_abs(Q + 0, Q + 0);
            fmpz_add_ui(Q + 0, Q + 0, 1);
            Rlen = 1 + n_randint(state, 2);
            for (i = 0; i < Rlen; i++)
                fmpz_randtest(R + i, state, 10);
        }

        for (i = 0; i < Plen; i++)
            _hs_set_int(Pi + i, P + i);
        for (i = 0; i < Qlen; i++)
            _hs_set_int(Qi + i, Q + i);
        for (i = 0; i < Rlen; i++)
            _hs_set_int(Ri + i, R + i);

        mp_real_init(T);
        mp_real_init(D);
        N2 = mp_real_hypgeom_sum(T, D, Pi, Plen, Qi, Qlen, Ri, Rlen, N, exact, n);

        if (N2 < N || (exact && N2 != N) || (N > 0 && 12 * (N2 - N) >= N))
        {
            flint_printf("FAIL: N = %wd, N' = %wd, exact = %d\n", N, N2, exact);
            flint_abort();
        }

        /* the exact sum of N2 terms as U / V, V = prod_(j<=k) Q(j), with
           Ar = prod_(j<k) R(j) (no gcds), and the largest term within a
           factor two, as 2^Mexp */
        fmpz_init(v);
        fmpz_init(w);
        fmpz_init(U);
        fmpz_init(V);
        fmpz_init(Ar);
        fmpz_one(V);
        fmpz_one(Ar);
        Mexp = WORD_MIN;
        for (k = 1; k <= N2; k++)
        {
            _hs_eval(v, P, Plen, k);
            _hs_eval(w, Q, Qlen, k);
            fmpz_mul(U, U, w);
            fmpz_mul(V, V, w);
            fmpz_mul(v, v, Ar);
            fmpz_add(U, U, v);
            if (!fmpz_is_zero(v))
                Mexp = FLINT_MAX(Mexp, (slong) fmpz_bits(v) - (slong) fmpz_bits(V) + 2);
            _hs_eval(v, R, Rlen, k);
            fmpz_mul(Ar, Ar, v);
        }
        fmpq_init(S);
        fmpq_set_fmpz_frac(S, U, V);

        arb_init(a);
        arb_init(b);
        if (N2 > 0)
        {
            mp_real_t q;
            mp_real_init(q);
            mp_real_div(q, T, D, n);
            mp_real_get_arb(a, q);
            mp_real_clear(q);
        }
        else
        {
            if (!mp_real_is_zero(T))
            {
                flint_printf("FAIL: empty sum\n");
                flint_abort();
            }
        }

        if (N2 > 0)
        {
            mag_t m;
            mag_init(m);
            if (Mexp != WORD_MIN)
                mag_set_ui_2exp_si(m, 1, Mexp - FLINT_BITS * (n - 1) + 24);
            if (!arb_contains_fmpq(a, S) || (Mexp != WORD_MIN && mag_cmp(arb_radref(a), m) > 0))
            {
                flint_printf("FAIL: N = %wd, N' = %wd, n = %wd, Plen %wd Qlen %wd Rlen %wd\n",
                    N, N2, n, Plen, Qlen, Rlen);
                flint_printf("a = "); arb_printd(a, 30); flint_printf("\n");
                arb_set_fmpq(b, S, FLINT_BITS * n + 64);
                flint_printf("S = "); arb_printd(b, 30); flint_printf("\n");
                flint_abort();
            }
            mag_clear(m);
        }

        for (i = 0; i < Plen; i++)
            flint_free((void *) Pi[i].d);
        for (i = 0; i < Qlen; i++)
            flint_free((void *) Qi[i].d);
        for (i = 0; i < Rlen; i++)
            flint_free((void *) Ri[i].d);
        for (i = 0; i < 4; i++)
        {
            fmpz_clear(P + i);
            fmpz_clear(Q + i);
            fmpz_clear(R + i);
        }
        fmpz_clear(v);
        fmpz_clear(w);
        fmpz_clear(U);
        fmpz_clear(V);
        fmpz_clear(Ar);
        fmpq_clear(S);
        mp_real_clear(T);
        mp_real_clear(D);
        arb_clear(a);
        arb_clear(b);
    }

    TEST_FUNCTION_END(state);
}
