/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_elliptic.h"
#include "acb_modular.h"

TEST_FUNCTION_START(acb_elliptic_pi_inc, state)
{
    slong iter;

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t z1, z2, n, m, r1, r2, r3, t;
        arb_t pi;
        fmpz_t k;
        slong prec1, prec2;
        int times_pi;

        acb_init(z1);
        acb_init(z2);
        acb_init(n);
        acb_init(m);
        acb_init(r1);
        acb_init(r2);
        acb_init(r3);
        acb_init(t);
        arb_init(pi);
        fmpz_init(k);

        prec1 = 2 + n_randint(state, 100);
        prec2 = 2 + n_randint(state, 100);
        times_pi = n_randint(state, 2);

        acb_randtest(n, state, 1 + n_randint(state, 500), 1 + n_randint(state, 10));
        acb_randtest(z1, state, 1 + n_randint(state, 500), 1 + n_randint(state, 10));
        acb_randtest(m, state, 1 + n_randint(state, 500), 1 + n_randint(state, 10));
        fmpz_randtest(k, state, 1 + n_randint(state, 100));
        arb_const_pi(pi, FLINT_MAX(prec1, prec2));

        acb_one(t);
        if ((acb_is_real(n) && acb_is_real(z1) && acb_is_real(m) &&
            arb_le(acb_realref(n), acb_realref(t)) &&
            arb_le(acb_realref(m), acb_realref(t))) || n_randint(state, 10) == 0)
        {
            if (n_randint(state, 2))
                arb_set_d(acb_realref(z1), -4.5 + n_randint(state, 10));

            /* test Pi(n,z,m) = Pi(n, z+pi k, m) - 2 k Pi(n, m) */

            if (times_pi)
            {
                if (n_randint(state, 2))
                {
                    acb_mul_arb(t, z1, pi, prec2);
                    acb_elliptic_pi_inc(r1, n, t, m, 0, prec1);
                }
                else
                {
                    acb_elliptic_pi_inc(r1, n, z1, m, 1, prec1);
                }

                if (n_randint(state, 2))
                {
                    acb_add_fmpz(z2, z1, k, prec2);
                    acb_elliptic_pi_inc(r2, n, z2, m, 1, prec2);
                }
                else
                {
                    acb_mul_arb(z2, z1, pi, prec2);
                    arb_addmul_fmpz(acb_realref(z2), pi, k, prec2);
                    acb_elliptic_pi_inc(r2, n, z2, m, 0, prec2);
                }
            }
            else
            {
                acb_elliptic_pi_inc(r1, n, z1, m, 0, prec1);

                acb_set(z2, z1);
                arb_addmul_fmpz(acb_realref(z2), pi, k, prec2);
                acb_elliptic_pi_inc(r2, n, z2, m, 0, prec2);
            }

            acb_set(r3, r2);
            acb_elliptic_pi(t, n, m, prec2);
            acb_mul_2exp_si(t, t, 1);
            acb_submul_fmpz(r3, t, k, prec2);

            if (!acb_overlaps(r1, r3))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("times_pi = %d\n\n", times_pi);
                flint_printf("z1 = "); acb_printd(z1, 30); flint_printf("\n\n");
                flint_printf("z2 = "); acb_printd(z2, 30); flint_printf("\n\n");
                flint_printf("n = "); acb_printd(n, 30); flint_printf("\n\n");
                flint_printf("m = "); acb_printd(m, 30); flint_printf("\n\n");
                flint_printf("k = "); fmpz_print(k); flint_printf("\n\n");
                flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
                flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
                flint_printf("r3 = "); acb_printd(r3, 30); flint_printf("\n\n");
                flint_abort();
            }

            /* test http://functions.wolfram.com/EllipticIntegrals/EllipticPi3/03/01/01/0006/ */
            /* Pi(n,z,n) = (n sin(2z) / (2 sqrt(1-n sin(z)^2)) - E(z,n))/(n-1) */
            acb_elliptic_pi_inc(r1, n, z1, n, times_pi, prec1);

            if (times_pi)
                acb_sin_pi(z2, z1, prec1);
            else
                acb_sin(z2, z1, prec1);
            acb_mul(z2, z2, z2, prec1);
            acb_mul(z2, z2, n, prec1);
            acb_sub_ui(z2, z2, 1, prec1);
            acb_neg(z2, z2);
            acb_rsqrt(z2, z2, prec1);
            acb_mul_2exp_si(z2, z2, -1);
            acb_mul_2exp_si(r3, z1, 1);
            if (times_pi)
                acb_sin_pi(r3, r3, prec1);
            else
                acb_sin(r3, r3, prec1);
            acb_mul(z2, z2, r3, prec1);
            acb_mul(z2, z2, n, prec1);
            acb_elliptic_e_inc(r3, z1, n, times_pi, prec1);
            acb_sub(z2, z2, r3, prec1);
            acb_sub_ui(r3, n, 1, prec1);
            acb_div(r2, z2, r3, prec1);

            if (!acb_overlaps(r1, r2))
            {
                flint_printf("FAIL: overlap Pi(n,z,n)\n\n");
                flint_printf("times_pi = %d\n\n", times_pi);
                flint_printf("z1 = "); acb_printd(z1, 30); flint_printf("\n\n");
                flint_printf("n = "); acb_printd(n, 30); flint_printf("\n\n");
                flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
                flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(z1);
        acb_clear(z2);
        acb_clear(n);
        acb_clear(m);
        acb_clear(r1);
        acb_clear(r2);
        acb_clear(r3);
        acb_clear(t);
        arb_clear(pi);
        fmpz_clear(k);
    }

    TEST_FUNCTION_END(state);
}
