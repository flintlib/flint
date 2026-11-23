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
#include "acb_hypgeom.h"

/* 1F1, U and 2F1 with large parameters, where the algorithm and the
   working precision are chosen from estimates: consistency between
   precisions, contiguous relations, containment for inexact input, and
   accuracy for exact real input in the domain of the integrals. */

static void
set_random(acb_t x, flint_rand_t state, double size, int cplx, int positive)
{
    double re, im = 0.0;

    re = size * (0.05 + n_randint(state, 1000) / 1000.0);
    if (!positive && n_randint(state, 4) == 0)
        re = -re;
    re += n_randint(state, 100) / 128.0;
    if (cplx)
        im = size * (n_randint(state, 2001) / 1000.0 - 1.0) * (n_randint(state, 2) ? 1.0 : 0.1);

    acb_set_d_d(x, re, im);
}

static void
eval(acb_t res, int f, const acb_t a, const acb_t b, const acb_t c, const acb_t z, slong prec)
{
    if (f == 0)
        acb_hypgeom_m(res, a, b, z, 0, prec);
    else if (f == 1)
        acb_hypgeom_u(res, a, b, z, prec);
    else
        acb_hypgeom_2f1(res, a, b, c, z, 0, prec);
}

TEST_FUNCTION_START(acb_hypgeom_large_parameters, state)
{
    slong iter;

    for (iter = 0; iter < 400 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, z, r1, r2, r3, t, u, am, ap;
        slong prec1, prec2;
        int f, cplx, easy, i;
        double size;

        acb_init(a); acb_init(b); acb_init(c); acb_init(z);
        acb_init(r1); acb_init(r2); acb_init(r3); acb_init(t); acb_init(u);
        acb_init(am); acb_init(ap);

        f = n_randint(state, 3);
        cplx = n_randint(state, 3) == 0;
        easy = !cplx && n_randint(state, 2);
        size = (n_randint(state, 2) ? 30.0 : 300.0) * (1 + n_randint(state, 4));
        prec1 = 30 + n_randint(state, 200);
        if (n_randint(state, 10) == 0)
            prec1 += 800;
        prec2 = prec1 + 10 + n_randint(state, 100);

        set_random(a, state, size, cplx, easy);
        set_random(b, state, size, cplx, easy);
        set_random(c, state, size, cplx, easy);
        set_random(z, state, size, cplx, 0);

        if (easy)
        {
            /* a >= 1, b - a >= 1 (and for 2F1: b >= 1, c - b >= 1) */
            acb_add_ui(a, a, 1, 200);
            acb_add(b, b, a, 200);
            acb_add_ui(b, b, 1, 200);
            acb_add(c, c, b, 200);
            acb_add_ui(c, c, 1, 200);
            if (f == 1)
                acb_abs(acb_realref(z), z, 200);
        }

        if (f == 2)
        {
            /* in the unit disc, or to the left */
            if (n_randint(state, 2))
                acb_set_d_d(z, 0.9 * (n_randint(state, 2001) / 1000.0 - 1.0),
                    cplx ? 0.4 * (n_randint(state, 2001) / 1000.0 - 1.0) : 0.0);
            else
                acb_set_d_d(z, -(double) n_randint(state, 20) - 0.5,
                    cplx ? (double) n_randint(state, 5) : 0.0);
        }

        eval(r1, f, a, b, c, z, prec1);
        eval(r2, f, a, b, c, z, prec2);

        if (!acb_overlaps(r1, r2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("f = %d, prec1 = %wd, prec2 = %wd\n", f, prec1, prec2);
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n");
            flint_printf("c = "); acb_printd(c, 30); flint_printf("\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
            flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n");
            flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n");
            flint_abort();
        }

        /* nearly full accuracy where the integrals apply */
        if (easy && acb_rel_accuracy_bits(r1) < prec1 - 25)
        {
            flint_printf("FAIL: accuracy\n\n");
            flint_printf("f = %d, prec1 = %wd, accuracy = %wd\n", f, prec1, acb_rel_accuracy_bits(r1));
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n");
            flint_printf("c = "); acb_printd(c, 30); flint_printf("\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
            flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n");
            flint_abort();
        }

        /* contiguous relations in a:
           (b-a) M(a-1) + (2a-b+z) M(a) - a M(a+1) = 0
           U(a-1) + (b-2a-z) U(a) + a (a-b+1) U(a+1) = 0
           (c-a) F(a-1) + (2a-c+(b-a)z) F(a) + a (z-1) F(a+1) = 0 */
        acb_sub_ui(am, a, 1, prec2);
        acb_add_ui(ap, a, 1, prec2);
        eval(r1, f, am, b, c, z, prec2);
        eval(r3, f, ap, b, c, z, prec2);

        if (f == 0)
        {
            acb_sub(t, b, a, prec2);
            acb_mul(r1, r1, t, prec2);
            acb_mul_2exp_si(t, a, 1);
            acb_sub(t, t, b, prec2);
            acb_add(t, t, z, prec2);
            acb_addmul(r1, r2, t, prec2);
            acb_submul(r1, r3, a, prec2);
        }
        else if (f == 1)
        {
            acb_mul_2exp_si(t, a, 1);
            acb_sub(t, b, t, prec2);
            acb_sub(t, t, z, prec2);
            acb_addmul(r1, r2, t, prec2);
            acb_sub(t, a, b, prec2);
            acb_add_ui(t, t, 1, prec2);
            acb_mul(t, t, a, prec2);
            acb_addmul(r1, r3, t, prec2);
        }
        else
        {
            acb_sub(t, c, a, prec2);
            acb_mul(r1, r1, t, prec2);
            acb_sub(t, b, a, prec2);
            acb_mul(t, t, z, prec2);
            acb_addmul_ui(t, a, 2, prec2);
            acb_sub(t, t, c, prec2);
            acb_addmul(r1, r2, t, prec2);
            acb_sub_ui(t, z, 1, prec2);
            acb_mul(t, t, a, prec2);
            acb_addmul(r1, r3, t, prec2);
        }

        if (!acb_contains_zero(r1))
        {
            flint_printf("FAIL: contiguous relation\n\n");
            flint_printf("f = %d, prec2 = %wd\n", f, prec2);
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n");
            flint_printf("c = "); acb_printd(c, 30); flint_printf("\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n");
            flint_printf("r = "); acb_printd(r1, 30); flint_printf("\n");
            flint_abort();
        }

        /* inexact input: the value at a point of the ball is contained */
        for (i = 0; i < 2; i++)
        {
            acb_t ab, bb, zb;
            mag_t m;

            acb_init(ab); acb_init(bb); acb_init(zb);
            mag_init(m);

            acb_set(ab, a);
            acb_set(bb, b);
            acb_set(zb, z);

            /* (inexact z only, or all of them) */
            acb_get_mag(m, z);
            mag_mul_2exp_si(m, m, -(slong) (10 + n_randint(state, prec1)));
            arb_add_error_mag(acb_realref(zb), m);
            if (cplx)
                arb_add_error_mag(acb_imagref(zb), m);

            if (i == 1)
            {
                acb_get_mag(m, a);
                mag_mul_2exp_si(m, m, -(slong) (10 + n_randint(state, prec1)));
                arb_add_error_mag(acb_realref(ab), m);
                acb_get_mag(m, b);
                mag_mul_2exp_si(m, m, -(slong) (10 + n_randint(state, prec1)));
                arb_add_error_mag(acb_realref(bb), m);
            }

            eval(r1, f, ab, bb, c, zb, prec1);

            if (!acb_overlaps(r1, r2))
            {
                flint_printf("FAIL: containment (inexact input, %d)\n\n", i);
                flint_printf("f = %d, prec1 = %wd, prec2 = %wd\n", f, prec1, prec2);
                flint_printf("a = "); acb_printd(ab, 30); flint_printf("\n");
                flint_printf("b = "); acb_printd(bb, 30); flint_printf("\n");
                flint_printf("c = "); acb_printd(c, 30); flint_printf("\n");
                flint_printf("z = "); acb_printd(zb, 30); flint_printf("\n");
                flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n");
                flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n");
                flint_abort();
            }

            acb_clear(ab); acb_clear(bb); acb_clear(zb);
            mag_clear(m);
        }

        acb_clear(a); acb_clear(b); acb_clear(c); acb_clear(z);
        acb_clear(r1); acb_clear(r2); acb_clear(r3); acb_clear(t); acb_clear(u);
        acb_clear(am); acb_clear(ap);
    }

    TEST_FUNCTION_END(state);
}
