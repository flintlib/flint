/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "fixed.h"

/* random sizes covering balanced, unbalanced and large shapes */
static void
_random_sizes(flint_rand_t state, mp_size_t * an, mp_size_t * bn)
{
    switch (n_randint(state, 7))
    {
        case 6:
            /* above the Newton cutoff of the dispatcher */
            *bn = 1024 + n_randint(state, 300);
            *an = *bn + 1024 + n_randint(state, 1500);
            break;
        case 0:
            *bn = 1 + n_randint(state, 10);
            *an = *bn + n_randint(state, 20);
            break;
        case 1:
            *bn = 1 + n_randint(state, 40);
            *an = *bn + n_randint(state, 400);
            break;
        case 2:
            *bn = 1 + n_randint(state, 500);
            *an = *bn + n_randint(state, 1500);
            break;
        case 3:
            *bn = 1 + n_randint(state, 500);
            *an = *bn + n_randint(state, 500);
            break;
        case 4:
            *bn = 300 + n_randint(state, 200);
            *an = *bn + 1 + n_randint(state, 4);
            break;
        default:
            *bn = 3 + n_randint(state, 60);
            *an = *bn + n_randint(state, 1000);
            break;
    }
}

TEST_FUNCTION_START(flint_mpn_tdiv_qr, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        mp_ptr a, b, q, r, q2, r2, binv;
        mp_size_t an, bn, n;
        int alg;

        _random_sizes(state, &an, &bn);
        n = an - bn + 1;

        a = flint_malloc(an * sizeof(mp_limb_t));
        b = flint_malloc(bn * sizeof(mp_limb_t));
        q = flint_malloc(n * sizeof(mp_limb_t));
        r = flint_malloc(bn * sizeof(mp_limb_t));
        q2 = flint_malloc(n * sizeof(mp_limb_t));
        r2 = flint_malloc(bn * sizeof(mp_limb_t));

        flint_mpn_rrandom(a, state, an);
        flint_mpn_rrandom(b, state, bn);
        if (b[bn - 1] == 0)
            b[bn - 1] = 1;

        /* sometimes make the quotient nearly exact or a divides b */
        if (n_randint(state, 4) == 0)
        {
            mp_ptr t = flint_malloc((n + bn) * sizeof(mp_limb_t));
            if (n >= bn)
                flint_mpn_mul(t, a, n, b, bn);
            else
                flint_mpn_mul(t, b, bn, a, n);
            flint_mpn_copyi(a, t, an);
            if (n_randint(state, 2))
                mpn_sub_1(a, a, an, n_randint(state, 3));
            else
                mpn_add_1(a, a, an, n_randint(state, 3));
            if (a[an - 1] == 0 && an > bn)
                a[an - 1] = 1;
            flint_free(t);
        }

        mpn_tdiv_qr(q2, r2, 0, a, an, b, bn);

        alg = n_randint(state, 7);
        if (alg == 0)
            flint_mpn_tdiv_qr(q, r, a, an, b, bn);
        else if (alg == 1)
            _flint_mpn_tdiv_qr_newton(q, r, a, an, b, bn);
        else if (alg == 2 && an > 2 * bn && bn >= 3)
            _flint_mpn_tdiv_qr_unbalanced(q, r, a, an, b, bn);
        else if (alg == 3 && bn >= 3 && an >= n + 2)
        {
            binv = flint_malloc((n + 4) * sizeof(mp_limb_t));
            fixed_inv_newton(binv, b, bn, n + 2);
            _flint_mpn_tdiv_qr_preinv(q, r, a, an, b, bn, binv, n + 2);
            flint_free(binv);
        }
        else if (alg == 4)
        {
            flint_mpn_tdiv_q(q, a, an, b, bn);
            flint_mpn_tdiv_r(r, a, an, b, bn);
        }
        else if (alg == 6)
        {
            _flint_mpn_tdiv_qr_preinvn(q, r, a, an, b, bn);
        }
        else
        {
            alg = 5;
            _flint_mpn_tdiv_qr_newton(q, NULL, a, an, b, bn);
            flint_mpn_copyi(r, r2, bn);
        }

        if (mpn_cmp(q, q2, n) != 0 || mpn_cmp(r, r2, bn) != 0)
            TEST_FUNCTION_FAIL("alg = %d, an = %wd, bn = %wd\na = %{ulong*}\nb = %{ulong*}\nq = %{ulong*}\nq2 = %{ulong*}\nr = %{ulong*}\nr2 = %{ulong*}\n",
                alg, an, bn, a, an, b, bn, q, n, q2, n, r, bn, r2, bn);

        flint_free(a);
        flint_free(b);
        flint_free(q);
        flint_free(r);
        flint_free(q2);
        flint_free(r2);
    }

    TEST_FUNCTION_END(state);
}
