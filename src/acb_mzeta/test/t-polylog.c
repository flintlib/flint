/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_mzeta.h"

TEST_FUNCTION_START(acb_mzeta_polylog, state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        acb_t z1, z2, z3, t;
        acb_ptr x;
        ulong s[8], e[8];
        slong r, i, prec = 2 + n_randint(state, 300);
        ulong N = 1 + n_randint(state, 6);

        acb_init(z1);
        acb_init(z2);
        acb_init(z3);
        acb_init(t);
        x = _acb_vec_init(8);

        switch (n_randint(state, 3))
        {
            case 0:
                /* depth 1 against acb_polylog: x with |x| <= 1, x != 1,
                   |1 - 1/x| >= 0.75 or |x| <= 1/2 */
                s[0] = 1 + n_randint(state, 5);
                if (n_randint(state, 2))
                {
                    arb_set_si(acb_realref(x), (slong) n_randint(state, 33) - 16);
                    arb_set_si(acb_imagref(x), (slong) n_randint(state, 33) - 16);
                    acb_mul_2exp_si(x, x, -5);
                }
                else
                {
                    e[0] = n_randint(state, N);
                    if (s[0] == 1 && e[0] == 0)
                        e[0] = (N > 1) ? 1 : 0;
                    if (s[0] == 1 && N == 1)
                        s[0] = 2;
                    {
                        fmpq_t q;
                        fmpq_init(q);
                        fmpq_set_si(q, 2 * e[0], N);
                        arb_sin_cos_pi_fmpq(acb_imagref(x), acb_realref(x), q, prec + 30);
                        fmpq_clear(q);
                    }
                }
                acb_mzeta_polylog(z1, s, x, 1, prec);
                acb_polylog_si(z2, s[0], x, prec);
                break;

            case 1:
                /* stuffle Li_a(x) Li_b(y) = Li_{a,b}(x,y) + Li_{b,a}(y,x) + Li_{a+b}(xy)
                   at roots of unity of order N <= 6 */
                {
                    ulong a = 2 + n_randint(state, 2), b = 1 + n_randint(state, 3);
                    ulong ab[2], ba[2], ex, ey, e2[2], c;
                    ex = n_randint(state, N);
                    ey = n_randint(state, N);
                    if (b == 1 && ey == 0)
                        b = 2;
                    ab[0] = a; ab[1] = b;
                    ba[0] = b; ba[1] = a;
                    acb_mzeta_mzv_colored(z1, &a, &ex, 1, N, prec);
                    acb_mzeta_mzv_colored(t, &b, &ey, 1, N, prec);
                    acb_mul(z1, z1, t, prec);
                    e2[0] = ex; e2[1] = ey;
                    acb_mzeta_mzv_colored(z2, ab, e2, 2, N, prec);
                    e2[0] = ey; e2[1] = ex;
                    acb_mzeta_mzv_colored(t, ba, e2, 2, N, prec);
                    acb_add(z2, z2, t, prec);
                    c = a + b;
                    e2[0] = (ex + ey) % N;
                    acb_mzeta_mzv_colored(t, &c, e2, 1, N, prec);
                    acb_add(z2, z2, t, prec);
                }
                break;

            default:
                /* MZVs as polylogarithms at 1 */
                {
                    arb_t m;
                    arb_init(m);
                    r = 1 + n_randint(state, 4);
                    for (i = 0; i < r; i++)
                    {
                        s[i] = 1 + n_randint(state, 3);
                        acb_one(x + i);
                        e[i] = 0;
                    }
                    if (s[0] < 2)
                        s[0] = 2;
                    acb_mzeta_polylog(z1, s, x, r, prec);
                    acb_mzeta_mzv(m, s, r, prec);
                    acb_set_arb(z2, m);
                    acb_mzeta_mzv_colored(z3, s, e, r, 1, prec);
                    if (!acb_overlaps(z1, z3))
                    {
                        flint_printf("FAIL (colored N = 1)\n");
                        flint_abort();
                    }
                    arb_clear(m);
                }
                break;
        }

        if (!acb_overlaps(z1, z2))
        {
            flint_printf("FAIL (iter %wd)\n", iter);
            flint_printf("prec = %wd, N = %wu\n", prec, N);
            acb_printd(z1, 30); flint_printf("\n");
            acb_printd(z2, 30); flint_printf("\n");
            flint_abort();
        }

        /* divergent: s_1 = 1, x_1 = 1 */
        s[0] = 1;
        s[1] = 2;
        acb_one(x);
        acb_one(x + 1);
        acb_mzeta_polylog(z1, s, x, 2, prec);
        if (acb_is_finite(z1))
        {
            flint_printf("FAIL (divergent)\n");
            flint_abort();
        }

        acb_clear(z1);
        acb_clear(z2);
        acb_clear(z3);
        acb_clear(t);
        _acb_vec_clear(x, 8);
    }

    TEST_FUNCTION_END(state);
}
