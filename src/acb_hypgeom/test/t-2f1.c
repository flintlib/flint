/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "acb.h"
#include "acb_poly.h"
#include "acb_hypgeom.h"
#include "acb_hypgeom/impl.h"

TEST_FUNCTION_START(acb_hypgeom_2f1, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, z, w1, w2, t;
        slong prec1, prec2;
        int reg1, reg2, ebits;
        int alg1, alg2;
        int flags1, flags2;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(z);
        acb_init(w1);
        acb_init(w2);
        acb_init(t);

        prec1 = 2 + n_randint(state, 300);
        prec2 = 2 + n_randint(state, 300);

        if (n_randint(state, 20) == 0)
            ebits = 30;
        else
            ebits = 5;

        flags1 = flags2 = 0;

        switch (n_randint(state, 3))
        {
            case 0:
                acb_set_si(a, n_randint(state, 500));
                acb_set_si(b, n_randint(state, 500));
                acb_set_si(c, n_randint(state, 10));
                break;
            case 1:
                acb_set_si(a, n_randint(state, 200) - 100);
                acb_set_si(b, n_randint(state, 200) - 100);
                acb_set_si(c, n_randint(state, 200) - 100);
                break;
            default:
                acb_randtest_param(a, state, 1 + n_randint(state, 400), 1 + n_randint(state, ebits));
                acb_randtest_param(b, state, 1 + n_randint(state, 400), 1 + n_randint(state, ebits));
                acb_randtest_param(c, state, 1 + n_randint(state, 400), 1 + n_randint(state, ebits));
        }

        /* TODO: test other cases (a-c, b-c), once implemented. */
        if (n_randint(state, 10) == 0)
        {
            acb_add_si(a, b, n_randint(state, 20) - 10, prec1);
            flags1 |= ACB_HYPGEOM_2F1_AB;
            if (n_randint(state, 2))
                flags2 |= ACB_HYPGEOM_2F1_AB;
        }

        if (n_randint(state, 10) == 0)
        {
            acb_add(c, a, b, prec1);
            acb_add_si(c, c, n_randint(state, 20) - 10, prec1);
            flags1 |= ACB_HYPGEOM_2F1_ABC;
            if (n_randint(state, 2))
                flags2 |= ACB_HYPGEOM_2F1_ABC;
        }

        acb_randtest_param(z, state, 1 + n_randint(state, 400), 1 + n_randint(state, ebits));
        acb_randtest(w1, state, 1 + n_randint(state, 400), 1 + n_randint(state, ebits));
        acb_randtest(w2, state, 1 + n_randint(state, 400), 1 + n_randint(state, ebits));

        if (n_randint(state, 2))
            arb_zero(acb_imagref(z));

        reg1 = n_randint(state, 2);
        reg2 = n_randint(state, 2);

        alg1 = n_randint(state, 10);
        alg2 = n_randint(state, 10);

        if (reg1) flags1 |= ACB_HYPGEOM_2F1_REGULARIZED;
        if (reg2) flags2 |= ACB_HYPGEOM_2F1_REGULARIZED;

        switch (alg1)
        {
            case 0:
                acb_hypgeom_2f1_direct(w1, a, b, c, z, reg1, prec1);
                break;
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
                acb_hypgeom_2f1_transform(w1, a, b, c, z, flags1, alg1, prec1);
                break;
            case 6:
                acb_hypgeom_2f1_corner(w1, a, b, c, z, reg1, prec1);
                break;
            default:
                acb_hypgeom_2f1(w1, a, b, c, z, flags1, prec1);
        }

        switch (alg2)
        {
            case 0:
                acb_hypgeom_2f1_direct(w2, a, b, c, z, reg2, prec2);
                break;
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
                acb_hypgeom_2f1_transform(w2, a, b, c, z, flags2, alg2, prec2);
                break;
            case 6:
                acb_hypgeom_2f1_corner(w2, a, b, c, z, reg2, prec2);
                break;
            default:
                acb_hypgeom_2f1(w2, a, b, c, z, flags2, prec2);
        }

        if (reg1 != reg2)
        {
            acb_rgamma(t, c, prec2);

            if (reg1)
                acb_mul(w2, w2, t, prec2);
            else
                acb_mul(w1, w1, t, prec2);
        }

        if (!acb_overlaps(w1, w2))
        {
            flint_printf("FAIL: consistency\n\n");
            flint_printf("iter = %wd, prec1 = %wd, prec2 = %wd\n\n", iter, prec1, prec2);
            flint_printf("alg1 = %d, alg2 = %d\n\n", alg1, alg2);
            flint_printf("reg1 = %d, reg2 = %d\n\n", reg1, reg2);
            flint_printf("flags1 = %d, flags2 = %d\n\n", flags1, flags2);
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_printf("w2 = "); acb_printd(w2, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(z);
        acb_clear(w1);
        acb_clear(w2);
        acb_clear(t);
    }

    /* A polynomial with a pole of the series beyond its last term:
       2F1(-100, 1/2; -200; 1/3).  All the terms are summed, as the tail
       bound of a truncated sum is infinite. */
    {
        acb_t a, b, c, z, w1, w2;
        slong prec;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(z);
        acb_init(w1);
        acb_init(w2);

        acb_set_si(a, -100);
        acb_set_d(b, 0.5);
        acb_set_si(c, -200);

        acb_set_ui(z, 1);
        acb_div_ui(z, z, 3, 1000);
        acb_hypgeom_2f1(w2, a, b, c, z, 0, 1000);

        for (prec = 128; prec <= 512; prec *= 2)
        {
            acb_set_ui(z, 1);
            acb_div_ui(z, z, 3, prec);
            acb_hypgeom_2f1(w1, a, b, c, z, 0, prec);

            if (!acb_overlaps(w1, w2) || acb_rel_accuracy_bits(w1) < prec - 10)
            {
                flint_printf("FAIL: polynomial with a pole beyond the last term\n\n");
                flint_printf("prec = %wd\n\n", prec);
                flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
                flint_printf("w2 = "); acb_printd(w2, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(z);
        acb_clear(w1);
        acb_clear(w2);
    }

    /* Parameters a, b > c > 0 with c - a and c - b negative integers:
       a polynomial in z / (z - 1) (times a power of 1 - z), whose terms
       have equal signs for 0 < z < 1 and cancel for z < 0.  The
       function is accurate at the precision of a double with exact
       input (more precision for the polynomial, or the integrals), and
       with input that is accurate to 53 bits (the integrals). */
    {
        acb_t a, b, c, z, w1, w2;
        slong i, acc;
        static const double zs[6] = { 0.5, 0.75, -1.0, -0.5, -10.0, 2.0 };

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(z);
        acb_init(w1);
        acb_init(w2);

        for (i = 0; i < 12; i++)
        {
            slong p = (i >= 10) ? 333 : 64;

            acb_set_ui(a, 423);
            acb_set_ui(b, 412);
            acb_set_ui(c, 153);
            acb_set_d(z, zs[i % 6]);

            if (i >= 6)
                acb_set_d(z, zs[2 + i % 2]);

            acb_hypgeom_2f1(w2, a, b, c, z, 0, 2500);

            if (i >= 8)
            {
                /* radii 2^-1000, z = -1 and -1/2: input which is more
                   accurate than the precision */
                mag_set_ui_2exp_si(arb_radref(acb_realref(a)), 1, -1000);
                mag_set_ui_2exp_si(arb_radref(acb_realref(b)), 1, -1000);
                mag_set_ui_2exp_si(arb_radref(acb_realref(c)), 1, -1000);
                mag_set_ui_2exp_si(arb_radref(acb_realref(z)), 1, -1000);
            }
            else if (i >= 6)
            {
                /* relative radii 2^-53, z = -1 and -1/2 */
                mag_set_ui_2exp_si(arb_radref(acb_realref(a)), 423, -53);
                mag_set_ui_2exp_si(arb_radref(acb_realref(b)), 412, -53);
                mag_set_ui_2exp_si(arb_radref(acb_realref(c)), 153, -53);
                mag_set_ui_2exp_si(arb_radref(acb_realref(z)), 1, -53);
            }

            acb_hypgeom_2f1(w1, a, b, c, z, 0, p);
            acc = acb_rel_accuracy_bits(w1);

            if (!acb_overlaps(w1, w2) || acb_rel_accuracy_bits(w2) < 2000
                || acc < ((i >= 8) ? p - 14 : (i >= 6) ? 30 : 50))
            {
                flint_printf("FAIL: a, b > c > 0, case %wd\n\n", i);
                flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
                flint_printf("w2 = "); acb_printd(w2, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(z);
        acb_clear(w1);
        acb_clear(w2);
    }

    /* The series of the limit cases of the connection formulas (power
       series in a parameter), and a series which is truncated before a
       pole of its lower parameter. */
    {
        acb_t a, b, c, z, w, w1;
        slong i, acc;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(z);
        acb_init(w);
        acb_init(w1);

        /* The terms of 2F1(2.03125, -1.203125; -81.8984375; 0.42...)
           fall below 2^-64 after 30 terms and stay there up to the
           pole at k = 82, but the tail bound at n = 30 is 2^52 times
           the term there: the series is summed past the pole. */
        acb_set_d(a, 2.03125);
        acb_set_d(b, -1.203125);
        acb_set_d(c, -81.8984375);
        acb_set_d(z, 0.42013549804687500);
        acb_hypgeom_2f1_direct(w1, a, b, c, z, 0, 1000);

        for (i = 0; i < 2; i++)
        {
            slong p = i ? 96 : 64;

            acb_hypgeom_2f1_direct(w, a, b, c, z, 0, p);
            acc = acb_rel_accuracy_bits(w);

            if (acb_rel_accuracy_bits(w1) < 900 || !acb_overlaps(w, w1) || acc < p - 12)
            {
                flint_printf("FAIL: truncation before a pole, prec %wd\n\n", p);
                flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
                flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        /* Limit cases of the connection formulas with complex numbers
           and |w| = 0.8 to 0.9: more terms than bits.  The forward
           recurrence for the power series lost all of the precision
           below 900 bits; rectangular or binary splitting does not. */
        {
            static const double lim[4][9] = {
                { 1.5, 0, 2.5, 0, 0.75, 0.5, 0.875, 0.75, 2 },
                { 0.25, 1, 3.25, 1, 1.5, 0, 0.25, 1, 3 },
                { 0.5, 0.25, 1.25, 0, 3.75, 0.25, 0.375, 0.625, 4 },
                { 2, 0, 0.5, -0.5, 2.5, -0.5, 0.75, 1, 5 },
            };

            for (i = 0; i < 4; i++)
            {
                slong p;

                acb_set_d_d(a, lim[i][0], lim[i][1]);
                acb_set_d_d(b, lim[i][2], lim[i][3]);
                acb_set_d_d(c, lim[i][4], lim[i][5]);
                acb_set_d_d(z, lim[i][6], lim[i][7]);

                acb_hypgeom_2f1(w1, a, b, c, z, 0, 2000);

                for (p = 64; p <= 512; p *= 2)
                {
                    acb_hypgeom_2f1_transform(w, a, b, c, z, 0, (int) lim[i][8], p);
                    acc = acb_rel_accuracy_bits(w);

                    if (acb_rel_accuracy_bits(w1) < 1980 || !acb_overlaps(w, w1) || acc < p - 30)
                    {
                        flint_printf("FAIL: limit case %wd, prec %wd\n\n", i, p);
                        flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
                        flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
                        flint_abort();
                    }
                }
            }
        }

        /* The number of terms of a power series with a parameter
           -m + x: the terms beyond m matter (they have the factor x
           in place of zero) and decrease like the others.  For
           2F1(1 + x, -999998 + x; 1 + x; -10^-6), 71 terms at 333
           bits (and not the largest number that is allowed, 3380). */
        {
            acb_poly_struct aa[2], bb[2];
            acb_poly_t zz;
            slong n;

            acb_poly_init(aa);
            acb_poly_init(aa + 1);
            acb_poly_init(bb);
            acb_poly_init(bb + 1);
            acb_poly_init(zz);

            acb_poly_set_coeff_si(aa, 0, 1);
            acb_poly_set_coeff_si(aa, 1, 1);
            acb_poly_set_coeff_si(aa + 1, 0, -999998);
            acb_poly_set_coeff_si(aa + 1, 1, 1);
            acb_poly_set_coeff_si(bb, 0, 1);
            acb_poly_set_coeff_si(bb, 1, 1);
            acb_poly_set_coeff_si(bb + 1, 0, 1);
            acb_set_d(z, -1e-6);
            acb_poly_set_acb(zz, z);

            n = acb_hypgeom_pfq_series_choose_n(aa, 2, bb, 2, zz, 2, 333);

            if (n < 40 || n > 200)
            {
                flint_printf("FAIL: number of terms with a parameter -m + x: %wd\n\n", n);
                flint_abort();
            }

            acb_poly_clear(aa);
            acb_poly_clear(aa + 1);
            acb_poly_clear(bb);
            acb_poly_clear(bb + 1);
            acb_poly_clear(zz);
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(z);
        acb_clear(w);
        acb_clear(w1);
    }

    /* The choice among the series, the connection formulas and the
       corner algorithm by estimates of their costs: fixed cases in
       which the algorithm chosen by z alone is slow or not accurate
       (the examples of issues 2964 and 2967; limit cases; z near 1
       with parameters for which the series at z has few terms; a
       series with a pole of the lower parameter beyond its terms that
       matter).  The results at 64 and 333 bits are accurate and
       consistent with the function at 2000 bits and with one formula
       (0: the series at z, 1 to 5 as in acb_hypgeom_2f1_transform, 6:
       the corner algorithm) at 3000 bits. */
    {
        static const double cases[14][9] = {
            /* a, b, c (real and imaginary parts), z, formula */
            { 5, 0, 7, 0, 33, 0, 0, 0.99, 1 },
            { 37, 0, 12, 0, 48, 0, 0.7523359028935359, -0.658779697028808, 4 },
            { 32, 0, 14, 0, 48, 0, 0.7298184519936511, 0.6836410074955939, 4 },
            { 423, 0, 412, 0, 153, 0, -23.23937447000604, 13.426370987830149, 2 },
            { 1, 0, 1000, 0, 2000, 0, 0.99, 0, 0 },
            { 5, 0, 5, 0, 12, 0, 0.711893, 0.340784, 0 },
            { 0.3125, 0, 0.6875, 0, 1.0000000000009094947, 0, 0.5, 1.0, 2 },
            { 1, 0, 1, 0, 1e6, 0, -1e6, 0, 3 },
            { 0.5, 100, 0.5, -100, 1, 0, -3, 0, 3 },
            { 2.03125, 0, -1.203125, 0, -81.8984375, 0, 0.42013549804687500, 0, 0 },
            { -7.84375, 0, 0.171875, 0, 15.2265625, 0, 0.980569, -0.011344, 4 },
            { 0.5, 0.25, 1.25, 0, 3.75, 0.25, 0.375, 0.625, 0 },
            { 26.53125, 0, -16.703125, 0, 7.7265625, 0, 22.314161, 4.562584, 2 },
            { 0.3125, 0, 0.25, 0, 0.1875, 0, 0.5, -0.875, 6 },
        };
        acb_t a, b, c, z, w, w1, w2;
        slong i, j, acc;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(z);
        acb_init(w);
        acb_init(w1);
        acb_init(w2);

        for (i = 0; i < 14; i++)
        {
            int alg = (int) cases[i][8];

            acb_set_d_d(a, cases[i][0], cases[i][1]);
            acb_set_d_d(b, cases[i][2], cases[i][3]);
            acb_set_d_d(c, cases[i][4], cases[i][5]);
            acb_set_d_d(z, cases[i][6], cases[i][7]);

            acb_hypgeom_2f1(w1, a, b, c, z, 0, 2000);

            if (alg == 0)
                acb_hypgeom_2f1_direct(w2, a, b, c, z, 0, 3000);
            else if (alg == 6)
                acb_hypgeom_2f1_corner(w2, a, b, c, z, 0, 3000);
            else
                acb_hypgeom_2f1_transform(w2, a, b, c, z, 0, alg, 3000);

            if (acb_rel_accuracy_bits(w1) < 1980 || acb_rel_accuracy_bits(w2) < 400
                || !acb_overlaps(w1, w2))
            {
                flint_printf("FAIL: choice of the algorithm, case %wd (references)\n\n", i);
                flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
                flint_printf("w2 = "); acb_printd(w2, 30); flint_printf("\n\n");
                flint_abort();
            }

            for (j = 0; j < 2; j++)
            {
                slong p = j ? 333 : 64;

                acb_hypgeom_2f1(w, a, b, c, z, 0, p);
                acc = acb_rel_accuracy_bits(w);

                if (!acb_overlaps(w, w1) || !acb_overlaps(w, w2) || acc < p - 12)
                {
                    flint_printf("FAIL: choice of the algorithm, case %wd, prec %wd\n\n", i, p);
                    flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
                    flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
                    flint_abort();
                }
            }
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(z);
        acb_clear(w);
        acb_clear(w1);
        acb_clear(w2);
    }

    /* The integrals where Euler's integral over [0, 1] does not apply:
       functions of the type of the Jacobi polynomials (a < 0 < c < b)
       for 0 < z < 1 and around it, where the series cancel by more
       bits than the parameters are large, and z on the cut (1, inf),
       with complex and with real saddle points.  The results at 64
       and 333 bits are accurate and consistent with a series at a
       precision which makes up for its cancellation. */
    {
        static const double cases[7][7] = {
            /* a, b, c, z, formula, precision of the reference */
            { -10000, 10001, 1, 0.5, 0, 0, 22000 },
            { -10000.5, 10001.25, 1.75, 0.3125, 0, 0, 22000 },
            { -8000.5, 9000.25, 2.5, 0.40625, 0.09375, 0, 22000 },
            { -3000, 3000.5, 0.5, 0.75, 0, 0, 16000 },
            { 1000, 2000, 5000, 2.5, 0, 2, 12000 },
            { 1000.25, 2000.5, 5000.125, 1.203125, 0, 5, 12000 },
            { 100.25, 200.5, 500.125, 1.203125, 0, 5, 3000 },
        };
        acb_t a, b, c, z, w, w2;
        slong i, j, acc;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(z);
        acb_init(w);
        acb_init(w2);

        for (i = 0; i < 7; i++)
        {
            int alg = (int) cases[i][5];

            acb_set_d(a, cases[i][0]);
            acb_set_d(b, cases[i][1]);
            acb_set_d(c, cases[i][2]);
            acb_set_d_d(z, cases[i][3], cases[i][4]);

            if (alg == 0)
                acb_hypgeom_2f1_direct(w2, a, b, c, z, 0, (slong) cases[i][6]);
            else
                acb_hypgeom_2f1_transform(w2, a, b, c, z, 0, alg, (slong) cases[i][6]);

            if (acb_rel_accuracy_bits(w2) < 400)
            {
                flint_printf("FAIL: loop and cut integrals, case %wd (reference)\n\n", i);
                flint_printf("w2 = "); acb_printd(w2, 30); flint_printf("\n\n");
                flint_abort();
            }

            for (j = 0; j < 2; j++)
            {
                slong p = j ? 333 : 64;

                acb_hypgeom_2f1(w, a, b, c, z, 0, p);
                acc = acb_rel_accuracy_bits(w);

                if (!acb_overlaps(w, w2) || acc < p - 12)
                {
                    flint_printf("FAIL: loop and cut integrals, case %wd, prec %wd\n\n", i, p);
                    flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
                    flint_printf("w2 = "); acb_printd(w2, 30); flint_printf("\n\n");
                    flint_abort();
                }
            }
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(z);
        acb_clear(w);
        acb_clear(w2);
    }

    /* The estimates for the connection formulas.  The function is at
       most the sum of the absolute values of the terms of a formula:
       the largest term Lam times the number of terms of that size (as
       in acb_hypgeom_2f1).  The size of the function where the
       integrand of Euler's integral has its maximum at an end point. */
    {
        acb_t a, b, c, z, w;
        arb_t x;
        slong iter2, n = 0;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(z);
        acb_init(w);
        arb_init(x);

        for (iter2 = 0; iter2 < 400 * 0.1 * flint_test_multiplier(); iter2++)
        {
            _acb_hypgeom_cd A, B, C, Z;
            _acb_hypgeom_est_2f1_conn_struct E;
            double lf;
            int which, cplx = n_randint(state, 2);

            A.re = (n_randint(state, 641) - 320.0) / 16.0 + 1.0 / 32.0;
            B.re = (n_randint(state, 641) - 320.0) / 16.0 + 3.0 / 64.0;
            C.re = (n_randint(state, 641) - 120.0) / 16.0 + 5.0 / 128.0;
            A.im = cplx ? (n_randint(state, 161) - 80.0) / 16.0 : 0.0;
            B.im = 0.0;
            C.im = cplx ? (n_randint(state, 161) - 80.0) / 16.0 : 0.0;
            Z.re = (n_randint(state, 801) - 400.0) / 64.0;
            Z.im = (n_randint(state, 801) - 400.0) / 64.0;

            acb_set_d_d(a, A.re, A.im);
            acb_set_d_d(b, B.re, B.im);
            acb_set_d_d(c, C.re, C.im);
            acb_set_d_d(z, Z.re, Z.im);

            acb_hypgeom_2f1(w, a, b, c, z, 0, 128);

            if (acb_rel_accuracy_bits(w) < 64)
                continue;

            acb_abs(x, w, 64);
            lf = log(arf_get_d(arb_midref(x), ARF_RND_NEAR));

            for (which = 2; which <= 5; which++)
            {
                double width;

                if (!_acb_hypgeom_est_2f1_conn(&E, which, A, B, C, Z) || E.limit)
                    continue;

                width = 2.0 * (FLINT_MAX(E.kpk[0], E.kpk[1]) + 1.0 / (-E.lw) + 6.0);
                n++;

                if (!(E.lw < 0.0) || !(E.Lam + log(width) > lf - 0.1))
                {
                    flint_printf("FAIL: largest term of a connection formula\n\n");
                    flint_printf("which = %d, Lam = %g, width = %g, log |f| = %g\n\n", which, E.Lam, width, lf);
                    flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
                    flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
                    flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
                    flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                    flint_abort();
                }
            }
        }

        if (n < 10 * flint_test_multiplier())
        {
            flint_printf("FAIL: few estimates for connection formulas (%wd)\n\n", n);
            flint_abort();
        }

        {
            _acb_hypgeom_cd A, B, C, Z;
            double f1, f2, g1, g2;

            A.re = 1000.0; A.im = 0.0;
            B.re = 1.0; B.im = 0.0;
            C.re = 2000.0; C.im = 0.0;
            Z.re = 0.99; Z.im = 0.0;

            /* 2F1(1000, 1; 2000; 0.99) = 1.98115... */
            _acb_hypgeom_est_2f1_logabs(&f1, &f2, A, B, C, Z, 1e300);
            _acb_hypgeom_est_2f1_logabs(&g1, &g2, B, A, C, Z, 1e300);

            if (!(fabs(f1 - 0.6837) < 0.5) || !(fabs(g1 - 0.6837) < 0.5))
            {
                flint_printf("FAIL: size of the function from an end point: %g, %g\n\n", f1, g1);
                flint_abort();
            }
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(z);
        acb_clear(w);
        arb_clear(x);
    }

    /* A ball z close to 1, where the function depends on 1 - z (the
       result contains the values at points of the ball, and it is as
       accurate as 1 - z allows, up to some bits), and large
       parameters with a series at z whose terms decrease by a factor
       of four or more from the first one (nearly the full accuracy). */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, z, zp, w, r;
        arf_t x;
        slong prec, rprec, k, rr, j;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(z);
        acb_init(zp);
        acb_init(w);
        acb_init(r);
        arf_init(x);

        prec = 40 + n_randint(state, 300);

        if (n_randint(state, 2))
        {
            /* a, b, c: multiples of 1/64 between -2 and 4 which are not integers */
            acb_set_si(a, 64 * ((slong) n_randint(state, 6) - 2) + 1 + n_randint(state, 63));
            acb_set_si(b, 64 * ((slong) n_randint(state, 6) - 2) + 1 + n_randint(state, 63));
            acb_set_si(c, 64 * ((slong) n_randint(state, 6) - 2) + 1 + n_randint(state, 63));
            acb_mul_2exp_si(a, a, -6);
            acb_mul_2exp_si(b, b, -6);
            acb_mul_2exp_si(c, c, -6);

            /* z = 1 - 2^-k with the radius 2^-rr, rr > k + 8 */
            k = 4 + n_randint(state, 40);
            rr = k + 8 + n_randint(state, prec);
            acb_one(z);
            arf_set_si_2exp_si(x, -1, -k);
            arb_add_arf(acb_realref(z), acb_realref(z), x, ARF_PREC_EXACT);
            mag_set_ui_2exp_si(arb_radref(acb_realref(z)), 1, -rr);

            rprec = prec + 2 * k + 100;
            acb_hypgeom_2f1(w, a, b, c, z, 0, prec);

            for (j = -2; j <= 2; j++)
            {
                acb_get_mid(zp, z);
                arf_set_mag(x, arb_radref(acb_realref(z)));
                arf_mul_si(x, x, j, ARF_PREC_EXACT, ARF_RND_DOWN);
                arf_mul_2exp_si(x, x, -1);
                arf_add(arb_midref(acb_realref(zp)), arb_midref(acb_realref(zp)), x, ARF_PREC_EXACT, ARF_RND_DOWN);

                acb_hypgeom_2f1(r, a, b, c, zp, 0, rprec);

                if (!acb_overlaps(w, r) || acb_rel_accuracy_bits(r) < prec)
                {
                    flint_printf("FAIL: point of a ball close to 1\n\n");
                    flint_printf("prec = %wd\n\n", prec);
                    flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
                    flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
                    flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
                    flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                    flint_printf("zp = "); acb_printd(zp, 30); flint_printf("\n\n");
                    flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
                    flint_printf("r = "); acb_printd(r, 30); flint_printf("\n\n");
                    flint_abort();
                }
            }

            /* (1 - z has rr - k accurate bits) */
            if (acb_rel_accuracy_bits(w) < FLINT_MIN(prec, rr - k) - 30)
            {
                flint_printf("FAIL: accuracy close to 1\n\n");
                flint_printf("prec = %wd, k = %wd, rr = %wd\n\n", prec, k, rr);
                flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
                flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
                flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
                flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
                flint_abort();
            }
        }
        else
        {
            /* |a| <= 2000, |b| <= 4, |z| <= 3/4,
               |c| >= 8 (|a| + 1) (|b| + 1) with Re(c) >= 0 */
            double aa, bb;

            arb_set_si(acb_realref(a), (slong) n_randint(state, 2001) - 1000);
            arb_set_si(acb_imagref(a), (slong) n_randint(state, 2001) - 1000);
            arb_set_si(acb_realref(b), (slong) n_randint(state, 65) - 32);
            arb_set_si(acb_imagref(b), (slong) n_randint(state, 65) - 32);
            acb_mul_2exp_si(b, b, -4);
            acb_div_ui(b, b, 3, prec);
            acb_get_mid(b, b);

            aa = 1415.0 + 1.0;
            bb = 3.0 + 1.0;
            arb_set_si(acb_realref(c), (slong) (8.0 * aa * bb) * (slong) n_randint(state, 2));
            arb_set_si(acb_imagref(c), (slong) (8.0 * aa * bb) * (1 + (slong) n_randint(state, 2)));
            if (n_randint(state, 2))
                arb_swap(acb_realref(c), acb_imagref(c));
            if (n_randint(state, 2))
                arb_neg(acb_imagref(c), acb_imagref(c));
            arb_set_d(acb_realref(z), 0.5);
            arb_add(acb_realref(c), acb_realref(c), acb_realref(z), prec);

            arb_set_si(acb_realref(z), (slong) n_randint(state, 65) - 32);
            arb_set_si(acb_imagref(z), (slong) n_randint(state, 65) - 32);
            acb_mul_2exp_si(z, z, -6);

            rprec = 2 * prec + 200;
            acb_hypgeom_2f1(w, a, b, c, z, 0, prec);
            acb_hypgeom_2f1(r, a, b, c, z, 0, rprec);

            if (!acb_overlaps(w, r) || acb_rel_accuracy_bits(w) < prec - 16)
            {
                flint_printf("FAIL: series with quickly decreasing terms\n\n");
                flint_printf("prec = %wd\n\n", prec);
                flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
                flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
                flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
                flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
                flint_printf("r = "); acb_printd(r, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(z);
        acb_clear(zp);
        acb_clear(w);
        acb_clear(r);
        arf_clear(x);
    }

    TEST_FUNCTION_END(state);
}
