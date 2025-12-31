/*
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2025 Guillaume Moroz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"
#include <math.h>

//TODO: make function add radius

static void set_radius(acb_ptr roots, slong len, double correction)
{
    slong i;
    mag_t rad, mag;
    mag_init(rad);
    mag_init(mag);
    for(i=0; i<len; i++) {
        acb_get_mag(mag, roots+i);
        mag_set_d(rad, (len-1)*correction);
        mag_mul(rad, rad, mag);
        acb_add_error_mag(roots+i, rad);
    }
    mag_clear(rad);
    mag_clear(mag);
}

static void check_polynomial(const acb_poly_t A, acb_ptr roots, double correction, slong prec, acb_poly_t B, acb_poly_t C)
{
    acb_t t;
    slong i, deg, prec2;
    deg = A->length - 1;
    set_radius(roots, deg, correction);

    prec2 = 500;
    if (correction < ldexp(1,-prec/2))
    {
        acb_poly_fit_length(B, 1);
        acb_set(B->coeffs, A->coeffs + deg);
        _acb_poly_set_length(B, 1);

        for (i = 0; i < deg; i++)
        {
            acb_poly_fit_length(C, 2);
            acb_one(C->coeffs + 1);
            acb_neg(C->coeffs + 0, roots + i);
            _acb_poly_set_length(C, 2);
            acb_poly_mul(B, B, C, prec2);
        }

        if (!acb_poly_contains(B, A))
        {
            flint_printf("FAIL: product does not equal polynomial\n");
            acb_poly_printd(A, 15); flint_printf("\n\n");
            acb_poly_printd(B, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_init(t);
        for (i = 0; i < deg; i++)
        {
            acb_poly_evaluate(t, A, roots + i, prec2);
            if (!acb_contains_zero(t))
            {
                flint_printf("FAIL: poly(root) does not contain zero\n");
                acb_poly_printd(A, 15); flint_printf("\n\n");
                acb_printd(roots + i, 15); flint_printf("\n\n");
                acb_printd(t, 15); flint_printf("\n\n");
                flint_abort();
            }
        }
        acb_clear(t);
    } else {
        flint_printf("FAIL: target precision was not reached by a large margin (%.3e >> %.3e)\n", correction, ldexp(1,-prec/2));
            acb_poly_printd(A, 15); flint_printf("\n\n");
            _acb_vec_sort_pretty(roots, deg);
            for(slong i=0; i<deg; i++) {
                flint_printf("%{acb}\n",roots+i);
            }
            flint_printf("\n\n");
            slong isolated = acb_poly_find_roots(roots, A, NULL, 150, prec2);
            _acb_vec_sort_pretty(roots, deg);
            for(slong i=0; i<deg; i++) {
                flint_printf("%{acb}\n",roots+i);
            }
            flint_printf("%ld isolated roots\n", isolated);
            flint_abort();
    }
}

TEST_FUNCTION_START(acb_poly_find_roots_double, state)
{
    slong iter;

    acb_poly_t A, B, C;
    acb_ptr roots;
    slong i, deg, prec;
    double correction;

    acb_poly_init(A);
    acb_poly_init(B);
    acb_poly_init(C);

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        prec = 10 + n_randint(state, 400);

        do {
            acb_poly_randtest(A, state, 50 + n_randint(state, 15), prec, 5);
            for(i = 0; i < A->length; i++) {
                arb_get_mid_arb(acb_realref(A->coeffs + i), acb_realref(A->coeffs + i));
                arb_get_mid_arb(acb_imagref(A->coeffs + i), acb_imagref(A->coeffs + i));
            }
        } while (A->length == 0 || acb_contains_zero(A->coeffs+A->length-1));
        deg = A->length - 1;

        roots = _acb_vec_init(deg);

        correction = _acb_poly_find_roots_double(roots, A->coeffs, NULL, A->length, 150, prec);
        check_polynomial(A, roots, correction, FLINT_MIN(53, prec), B, C);

        _acb_vec_clear(roots, deg);
    }

    /* Check rescaling in an extreme case */
    slong hdeg = 100;

    acb_poly_fit_length(A, 2*hdeg+1);
    _acb_vec_zero(A->coeffs, 2*hdeg+1);
    acb_poly_set_coeff_si(A, 0, 1);
    acb_poly_set_coeff_si(A, 1, 1);
    acb_poly_set_coeff_si(A, hdeg, -1);
    acb_mul_2exp_si(A->coeffs+hdeg, A->coeffs+hdeg, 2000);
    acb_poly_set_coeff_si(A, 2*hdeg, 1);
    _acb_poly_set_length(A, 2*hdeg+1);

    deg = A->length - 1;

    roots = _acb_vec_init(deg);

    correction = _acb_poly_find_roots_double(roots, A->coeffs, NULL, A->length, 150, 53);
    check_polynomial(A, roots, correction, 53, B, C);

    _acb_vec_clear(roots, deg);

    /* Check specific polynomial with clustered roots */
    double f[51][2] = {
        {-1.0000000000000000, 0.99998474132734927},
        {-0.69891357421875000, 1.0000000000000000},
        {-9.5367431640625000e-7, 0.49975585960782712},
        {0.50000000000000000, 16776672.031249881},
        {0.49999999999999645, 0.0078125000056843410},
        {-0.54034423828125000, 0.24999999999999999},
        {127.99926757812318, 0.99999994039535523},
        {1.0000000000000000, 4.0000000000000000},
        {0.12499999906867743, 1.8596252012030163},
        {-1.9999999944102465, 0.99999999999909051},
        {0, 0.49951171875011368},
        {4.7683715820312479e-7, 65535.999999985099},
        {-0.50000000000000000, 1.9999389648437518},
        {0.57812500000000000, 0.00018310593440773459},
        {-24.000011444441043, 0.49999999998544809},
        {-1.8626451492309570e-9, 0.50000000000000000},
        {0.031250000000000000, 3.9999999998908322},
        {0.99999976158142090, 3.0517578123223643e-5},
        {20137.148437500000, 0.0078124999999998890},
        {2.0000000000000000, 2.9802324164049187e-8},
        {0.12499999999818101, 0.97653959045974403},
        {-0.12499999999998579, 0.50000000000000000},
        {-3.9687652587890625, 0.062499999985448085},
        {0.99999976158142179, 0.25000000000000000},
        {-3.8146972656250000e-6, 134217727.99993908},
        {-0.00048826634929355350, 1.9995117187572724},
        {15.750122069381632, 0.0078124701976776123},
        {-0.93923187255859375, 0.49999999627561919},
        {2.3841857909809306e-7, 5.6916985511779785},
        {-0.00047302246093750000, 0.46899414045037699},
        {0.24999999999317879, 4.7683715820312161e-7},
        {0.99999999999998623, 0.49999999999272404},
        {-8.0000000000000000, 0.0078124999999999445},
        {0.96875023841494134, 0.013870531693100929},
        {-0.99999998509883881, 1.7980164604199005},
        {0.52285766601562500, 0.00012207031250000000},
        {4095.9687500000000, 0.99993896927071590},
        {-0.0039024353045533609, 5.3200810256286913e-5},
        {-0.49999999813735485, 0.50000000000000000},
        {0.24951195716857888, 2.0000000000000000},
        {0.75000762939453125, 65535.999984741211},
        {0.00012207031250000000, 3.9999923706054829},
        {32768.000000000000, 0.50000000000000000},
        {-0.49999999994179323, 65536.000000000000},
        {-33553408.000000000, 7.9384536724573989},
        {-1.7573979979590426e-9, 1.9999999850997483},
        {-15.999999046325712, 0.031242370605468750},
        {-9.5320867998815484e-7, 1023.9999999995344},
        {134217727.99998474, 16.000000000000000},
        {8183.7500000000000, 128.00000000000000},
        {0.49999999999977263, 1.8568480868452752e-6}};

    acb_poly_fit_length(A, 51);
    for(i=0; i<51; i++) {
        acb_set_d_d(A->coeffs + i, f[i][0], f[i][1]);
    }
    _acb_poly_set_length(A, 51);

    deg = A->length-1;

    roots = _acb_vec_init(deg);

    correction = _acb_poly_find_roots_double(roots, A->coeffs, NULL, A->length, 150, 53);
    check_polynomial(A, roots, correction, 53, B, C);

    _acb_vec_clear(roots, deg);

    acb_poly_clear(A);
    acb_poly_clear(B);
    acb_poly_clear(C);

    TEST_FUNCTION_END(state);
}

TEST_FUNCTION_START(cd_poly_refine_roots, state)
{
    slong iter;

    mag_t rad;
    acb_poly_t A, B, C;
    acb_ptr roots;
    double *p_r, *p_i, *z_r, *z_i, *vp_r, *vp_i, *wdk_r, *wdk_i, *p, *z;
    slong i, deg, len, prec, max_deg;
    double correction;

    max_deg=100;
    mag_init(rad);
    acb_poly_init(A);
    acb_poly_init(B);
    acb_poly_init(C);
    roots = _acb_vec_init(max_deg);

    p_r = flint_malloc((max_deg+1)*sizeof(double));
    p_i = flint_malloc((max_deg+1)*sizeof(double));
    z_r = flint_malloc(max_deg*sizeof(double));
    z_i = flint_malloc(max_deg*sizeof(double));
    vp_r = flint_malloc(max_deg*sizeof(double));
    vp_i = flint_malloc(max_deg*sizeof(double));
    wdk_r = flint_malloc(max_deg*sizeof(double));
    wdk_i = flint_malloc(max_deg*sizeof(double));
    p = flint_malloc(2*(max_deg+1)*sizeof(double));
    z = flint_malloc(2*max_deg*sizeof(double));

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        prec = 10 + n_randint(state, 400);

        do {
            acb_poly_randtest(A, state, 50 + n_randint(state, 15), prec, 5);
            for(i = 0; i < A->length; i++) {
                p_r[i] = arf_get_d(arb_midref(acb_realref(A->coeffs + i)), ARF_RND_NEAR);
                p_i[i] = arf_get_d(arb_midref(acb_imagref(A->coeffs + i)), ARF_RND_NEAR);
                p[2*i]   = p_r[i];
                p[2*i+1] = p_i[i];
                arb_get_mid_arb(acb_realref(A->coeffs + i), acb_realref(A->coeffs + i));
                arb_get_mid_arb(acb_imagref(A->coeffs + i), acb_imagref(A->coeffs + i));
            }
        } while (A->length == 0 || acb_contains_zero(A->coeffs+A->length-1) || acb_contains_zero(A->coeffs));

        len = A->length;
        deg = A->length - 1;

        correction = 1;
        cd_poly_roots_initial_values(z_r, z_i, p_r, p_i, len, NULL, deg);
        for(i=0; i<deg; i++) {
            z[2*i]   = z_r[i];
            z[2*i+1] = z_i[i];
        }
        for(i=0; i<1000 && (deg-1)*correction > ldexp(1,-prec) && correction > ldexp(1,-50); i++) {
            correction = cd_poly_refine_roots_with_pivot(z, p, len, 0x1p-10);
        }
        for(i=0; i<10; i++) {
            correction = cd_poly_refine_roots(z, p, len, 0x1p-10);
        }
        for(i=0; i<deg; i++) {
            acb_set_d_d(roots + i, z[2*i], z[2*i+1]);
            mag_set_d(rad, (deg-1)*correction*hypot(z[2*i], z[2*i+1]));
            acb_add_error_mag(roots+i, rad);
        }

        /* Check that initial roots are correctly handled */
        correction = cd_poly_find_roots(z, p, z, len, 1, 0x1p-53);
        correction = _acb_poly_find_roots_double(roots, A->coeffs, roots, A->length, 1, prec);

        check_polynomial(A, roots, correction, FLINT_MIN(53, prec), B, C);
    }
    _acb_vec_clear(roots, max_deg);

    acb_poly_clear(A);
    acb_poly_clear(B);
    acb_poly_clear(C);
    mag_clear(rad);
    flint_free(p_r);
    flint_free(p_i);
    flint_free(z_r);
    flint_free(z_i);
    flint_free(vp_r);
    flint_free(vp_i);
    flint_free(wdk_r);
    flint_free(wdk_i);
    flint_free(p);
    flint_free(z);

    TEST_FUNCTION_END(state);
}

