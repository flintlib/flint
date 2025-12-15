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

void check_polynomial(const acb_poly_t A, acb_ptr roots, double correction, slong prec)
{
    acb_poly_t B;
    acb_poly_t C;
    acb_t t;
    slong i, deg, prec2;
    deg = A->length - 1;

    acb_init(t);
    acb_poly_init(B);
    acb_poly_init(C);
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

    acb_clear(t);
    acb_poly_clear(B);
    acb_poly_clear(C);
}

TEST_FUNCTION_START(acb_poly_find_roots_double, state)
{
    slong iter;

    acb_poly_t A;
    acb_ptr roots;
    slong i, deg, prec;
    double correction;

    acb_poly_init(A);

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
        check_polynomial(A, roots, correction, FLINT_MIN(53, prec));

        _acb_vec_clear(roots, deg);
    }

    /* Check rescaling in an extreme case */
    slong hdeg = 90;
    acb_poly_clear(A);
    acb_poly_init(A);
    acb_poly_fit_length(A, 2*hdeg+1);
    acb_poly_set_coeff_si(A, 0, 1);
    acb_poly_set_coeff_si(A, 1, 1);
    acb_poly_set_coeff_si(A, hdeg, -1);
    acb_mul_2exp_si(A->coeffs+hdeg, A->coeffs+hdeg, 1800);
    acb_poly_set_coeff_si(A, 2*hdeg, 1);
    _acb_poly_set_length(A, 2*hdeg+1);

    deg = A->length - 1;

    roots = _acb_vec_init(deg);

    correction = _acb_poly_find_roots_double(roots, A->coeffs, NULL, A->length, 150, 53);
    check_polynomial(A, roots, correction, 53);
    
    _acb_vec_clear(roots, deg);

    acb_poly_clear(A);

    TEST_FUNCTION_END(state);
}

TEST_FUNCTION_START(cd_poly_refine_roots, state)
{
    slong iter;

    mag_t rad;
    acb_poly_t A;
    acb_ptr roots;
    double *p_r, *p_i, *z_r, *z_i, *vp_r, *vp_i, *wdk_r, *wdk_i, *p, *z;
    slong i, deg, len, prec, max_deg;
    double correction;

    max_deg=100;
    mag_init(rad);
    acb_poly_init(A);
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
        for(i=0; i<10 && (deg-1)*correction > ldexp(1,-prec) && correction > ldexp(1,-50); i++) {
            correction = cd_poly_refine_roots(z, p, len, 0x1p-10);
        }
        for(i=0; i<deg; i++) {
            acb_set_d_d(roots + i, z[2*i], z[2*i+1]);
            mag_set_d(rad, (deg-1)*correction*hypot(z[2*i], z[2*i+1]));
            acb_add_error_mag(roots+i, rad);
        }

        /* Check initial roots are correctly handled */
        correction = cd_poly_find_roots(z, p, z, len, 1, 0x1p-53);
        correction = _acb_poly_find_roots_double(roots, A->coeffs, roots, A->length, 1, prec);

        check_polynomial(A, roots, correction, FLINT_MIN(53, prec));
    }
    _acb_vec_clear(roots, max_deg);
 
    acb_poly_clear(A);
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

