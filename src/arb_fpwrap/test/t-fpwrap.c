/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include <mpfr.h>
#include "double_extras.h"
#include "arb.h"
#include "arb_fpwrap.h"

#define CHECK_DOUBLE(fcall) \
    do { \
        int fail = fcall; \
        if (fail || res != res) \
        { \
            flint_printf("FAIL\n"); \
            flint_printf("%d\n", __LINE__); \
            flint_abort(); \
        } \
    } while (0)

#define CHECK_CDOUBLE(fcall) \
    do { \
        int fail = fcall; \
        if (fail || cres.real != cres.real || cres.imag != cres.imag) \
        { \
            flint_printf("FAIL\n"); \
            flint_printf("%d\n", __LINE__); \
            flint_abort(); \
        } \
    } while (0)

TEST_FUNCTION_START(arb_fpwrap, state)
{
    slong iter;

    /* correct rounding test */
    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        mpfr_t t;
        double x, y, z;

        mpfr_init2(t, 53);

        x = d_randtest(state) + n_randint(state, 100);

        mpfr_set_d(t, x, MPFR_RNDN);

        switch (n_randint(state, 4))
        {
            case 0:
                arb_fpwrap_double_log1p(&y, x, FPWRAP_CORRECT_ROUNDING);
                mpfr_log1p(t, t, MPFR_RNDN);
                break;
            case 1:
                arb_fpwrap_double_sqrt(&y, x, FPWRAP_CORRECT_ROUNDING);
                mpfr_sqrt(t, t, MPFR_RNDN);
                break;
            case 2:
                arb_fpwrap_double_exp(&y, x, FPWRAP_CORRECT_ROUNDING);
                mpfr_exp(t, t, MPFR_RNDN);
                break;
            default:
                arb_fpwrap_double_sin(&y, x, FPWRAP_CORRECT_ROUNDING);
                mpfr_sin(t, t, MPFR_RNDN);
                break;
        }

        z = mpfr_get_d(t, MPFR_RNDN);

        if (z != y)
        {
            flint_printf("FAIL: correct rounding\n\n");
            flint_abort();
        }

        mpfr_clear(t);
    }

    {
        double a[1], b[2];
        double z, y;

        z = 1.75;
        a[0] = 0.25;
        b[0] = 1.5;
        b[1] = -2.125;

        arb_fpwrap_double_hypgeom_pfq(&y, a, 1, b, 2, z, 0, 0);

        if (fabs(y - 0.68910385124070327187) > 1e-16)
        {
            flint_printf("FAIL: value 1\n\n");
            flint_abort();
        }

        arb_fpwrap_double_hypgeom_pfq(&y, a, 1, b, 2, z, 1, 0);

        if (fabs(y - (-0.21324224371323783595)) > 1e-16)
        {
            flint_printf("FAIL: value 2\n\n");
            flint_abort();
        }
    }

    {
        complex_double y, z;

        arb_fpwrap_cdouble_zeta_zero(&y, 1, FPWRAP_CORRECT_ROUNDING);
        arb_fpwrap_cdouble_zeta(&z, y, 0);

        if (fabs(z.real - (-1.0483650805588237388e-16)) > 1e-31)
        {
            flint_printf("FAIL: value 3\n\n");
            flint_abort();
        }

        if (fabs(z.imag - 6.5852592776051578103e-16) > 1e-31)
        {
            flint_printf("FAIL: value 4\n\n");
            flint_abort();
        }
    }

    {
        complex_double x, y;

        x.real = 1.0;
        x.imag = 1e-100;

        arb_fpwrap_cdouble_erf(&y, x, FPWRAP_ACCURATE_PARTS);

        if (fabs(y.imag - 4.1510749742059471164e-101) > 1e-116)
        {
            flint_printf("FAIL: value 5\n\n");
            flint_abort();
        }
    }

    {
        double y;

        arb_fpwrap_double_hypgeom_2f1(&y, 0.1, 0.2, 0.3, 0.4, 0, FPWRAP_CORRECT_ROUNDING);

        if (fabs(y - 1.0341794015503748492) > 1e-16)
        {
            flint_printf("FAIL: value 6\n\n");
            flint_abort();
        }

        arb_fpwrap_double_hypgeom_2f1(&y, 0.1, 0.2, 0.3, 0.4, 1, FPWRAP_CORRECT_ROUNDING);

        if (fabs(y - 0.34569799520143110351) > 1e-16)
        {
            flint_printf("FAIL: value 6\n\n");
            flint_abort();
        }
    }

    {
        double x, y, z;
        complex_double cx, cy, cz, ctau;
        double res, res2;
        complex_double cres;
        int flags;

        x = 0.25;
        y = 0.625;
        z = 0.75;

        cx.real = 0.25;
        cx.imag = 0.125;
        cy.real = 0.5;
        cy.imag = 0.625;
        cz.real = 0.75;
        cz.imag = 0.125;
        ctau.real = 0.25;
        ctau.imag = 1.0;

        flags = 0;

        CHECK_DOUBLE(arb_fpwrap_double_exp(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_exp(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_expm1(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_expm1(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_log(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_log(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_log1p(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_log1p(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_pow(&res, x, y, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_pow(&cres, cx, cy, flags));

        CHECK_DOUBLE(arb_fpwrap_double_sqrt(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_sqrt(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_rsqrt(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_rsqrt(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_cbrt(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_cbrt(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_sin(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_sin(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_cos(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_cos(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_tan(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_tan(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_cot(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_cot(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_sec(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_sec(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_csc(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_csc(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_sinc(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_sinc(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_sin_pi(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_sin_pi(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_cos_pi(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_cos_pi(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_tan_pi(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_tan_pi(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_cot_pi(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_cot_pi(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_sinc_pi(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_sinc_pi(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_asin(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_asin(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_acos(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_acos(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_atan(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_atan(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_atan2(&res, x, y, flags));

        CHECK_DOUBLE(arb_fpwrap_double_asinh(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_asinh(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_acosh(&res, 1.0 + x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_acosh(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_atanh(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_atanh(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_lambertw(&res, x, 0, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_lambertw(&cres, cx, 0, flags));

        CHECK_DOUBLE(arb_fpwrap_double_lambertw(&res, -0.2, -1, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_lambertw(&cres, cx, -1, flags));

        CHECK_DOUBLE(arb_fpwrap_double_rising(&res, x, y, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_rising(&cres, cx, cy, flags));

        CHECK_DOUBLE(arb_fpwrap_double_gamma(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_gamma(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_rgamma(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_rgamma(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_lgamma(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_lgamma(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_digamma(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_digamma(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_zeta(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_zeta(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_hurwitz_zeta(&res, x, z, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_hurwitz_zeta(&cres, cx, cz, flags));

        CHECK_DOUBLE(arb_fpwrap_double_lerch_phi(&res, x, y, z, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_lerch_phi(&cres, cx, cy, cz, flags));

        CHECK_DOUBLE(arb_fpwrap_double_barnes_g(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_barnes_g(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_log_barnes_g(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_log_barnes_g(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_polygamma(&res, x, z, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_polygamma(&cres, cx, cz, flags));

        CHECK_DOUBLE(arb_fpwrap_double_polylog(&res, x, z, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_polylog(&cres, cx, cz, flags));

        CHECK_CDOUBLE(arb_fpwrap_cdouble_dirichlet_eta(&cres, cz, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_riemann_xi(&cres, cz, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_hardy_theta(&cres, cz, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_hardy_z(&cres, cz, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_zeta_zero(&cres, 2, flags));

        CHECK_DOUBLE(arb_fpwrap_double_erf(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_erf(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_erfc(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_erfc(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_erfi(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_erfi(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_erfinv(&res, x, flags));
        CHECK_DOUBLE(arb_fpwrap_double_erfcinv(&res, x, flags));

        CHECK_DOUBLE(arb_fpwrap_double_fresnel_s(&res, x, 0, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_fresnel_s(&cres, cx, 0, flags));

        CHECK_DOUBLE(arb_fpwrap_double_fresnel_c(&res, x, 0, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_fresnel_c(&cres, cx, 0, flags));

        CHECK_DOUBLE(arb_fpwrap_double_gamma_upper(&res, x, z, 0, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_gamma_upper(&cres, cx, cz, 0, flags));

        CHECK_DOUBLE(arb_fpwrap_double_gamma_lower(&res, x, z, 0, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_gamma_lower(&cres, cx, cz, 0, flags));

        CHECK_DOUBLE(arb_fpwrap_double_beta_lower(&res, x, y, z, 0, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_beta_lower(&cres, cx, cy, cz, 0, flags));

        CHECK_DOUBLE(arb_fpwrap_double_exp_integral_e(&res, x, z, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_exp_integral_e(&cres, cx, cz, flags));

        CHECK_DOUBLE(arb_fpwrap_double_exp_integral_ei(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_exp_integral_ei(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_sin_integral(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_sin_integral(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_cos_integral(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_cos_integral(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_sinh_integral(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_sinh_integral(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_cosh_integral(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_cosh_integral(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_log_integral(&res, x, 0, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_log_integral(&cres, cx, 0, flags));

        CHECK_DOUBLE(arb_fpwrap_double_dilog(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_dilog(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_bessel_j(&res, y, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_bessel_j(&cres, cy, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_bessel_y(&res, y, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_bessel_y(&cres, cy, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_bessel_i(&res, y, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_bessel_i(&cres, cy, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_bessel_k(&res, y, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_bessel_k(&cres, cy, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_bessel_k_scaled(&res, y, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_bessel_k_scaled(&cres, cy, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_airy_ai(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_airy_ai(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_airy_ai_prime(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_airy_ai_prime(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_airy_bi(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_airy_bi(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_airy_bi_prime(&res, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_airy_bi_prime(&cres, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_airy_ai_zero(&res, 1, flags));
        CHECK_DOUBLE(arb_fpwrap_double_airy_ai_prime_zero(&res, 1, flags));
        CHECK_DOUBLE(arb_fpwrap_double_airy_bi_zero(&res, 1, flags));
        CHECK_DOUBLE(arb_fpwrap_double_airy_bi_prime_zero(&res, 1, flags));

        CHECK_DOUBLE(arb_fpwrap_double_coulomb_f(&res, y, z, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_coulomb_f(&cres, cy, cz, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_coulomb_g(&res, y, z, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_coulomb_g(&cres, cy, cz, cx, flags));

        CHECK_CDOUBLE(arb_fpwrap_cdouble_coulomb_hpos(&cres, cy, cz, cx, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_coulomb_hneg(&cres, cy, cz, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_chebyshev_t(&res, y, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_chebyshev_t(&cres, cy, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_chebyshev_u(&res, y, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_chebyshev_u(&cres, cy, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_jacobi_p(&res, y, x, y, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_jacobi_p(&cres, cy, cx, cy, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_gegenbauer_c(&res, y, z, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_gegenbauer_c(&cres, cy, cz, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_laguerre_l(&res, y, z, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_laguerre_l(&cres, cy, cz, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_hermite_h(&res, y, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_hermite_h(&cres, cy, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_legendre_p(&res, y, z, x, 0, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_legendre_p(&cres, cy, cz, cx, 0, flags));

        CHECK_DOUBLE(arb_fpwrap_double_legendre_q(&res, y, z, x, 0, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_legendre_q(&cres, cy, cz, cx, 0, flags));

        CHECK_DOUBLE(arb_fpwrap_double_legendre_root(&res, &res2, 2, 1, flags));

        CHECK_CDOUBLE(arb_fpwrap_cdouble_spherical_y(&cres, 2, 1, cx, cy, flags));

        CHECK_DOUBLE(arb_fpwrap_double_hypgeom_0f1(&res, y, x, 0, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_hypgeom_0f1(&cres, cy, cx, 0, flags));

        CHECK_DOUBLE(arb_fpwrap_double_hypgeom_1f1(&res, x, y, x, 0, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_hypgeom_1f1(&cres, cx, cy, cx, 0, flags));

        CHECK_DOUBLE(arb_fpwrap_double_hypgeom_u(&res, x, y, x, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_hypgeom_u(&cres, cx, cy, cx, flags));

        CHECK_DOUBLE(arb_fpwrap_double_hypgeom_2f1(&res, x, y, z, x, 0, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_hypgeom_2f1(&cres, cx, cy, cz, cx, 0, flags));

        CHECK_DOUBLE(arb_fpwrap_double_hypgeom_pfq(&res, &x, 1, &y, 1, z, 0, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_hypgeom_pfq(&cres, &cx, 1, &cy, 1, cz, 0, flags));

        CHECK_DOUBLE(arb_fpwrap_double_agm(&res, x, y, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_agm(&cres, cx, cy, flags));

        CHECK_CDOUBLE(arb_fpwrap_cdouble_elliptic_k(&cres, cz, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_elliptic_e(&cres, cz, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_elliptic_pi(&cres, cy, cz, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_elliptic_f(&cres, cx, cz, 0, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_elliptic_e_inc(&cres, cx, cz, 0, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_elliptic_pi_inc(&cres, cx, cy, cz, 0, flags));

        CHECK_CDOUBLE(arb_fpwrap_cdouble_elliptic_rf(&cres, cx, cy, cz, 0, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_elliptic_rg(&cres, cx, cy, cz, 0, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_elliptic_rj(&cres, cx, cy, cz, cx, 0, flags));

        CHECK_CDOUBLE(arb_fpwrap_cdouble_elliptic_p(&cres, cz, ctau, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_elliptic_p_prime(&cres, cz, ctau, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_elliptic_inv_p(&cres, cz, ctau, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_elliptic_zeta(&cres, cz, ctau, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_elliptic_sigma(&cres, cz, ctau, flags));

        CHECK_CDOUBLE(arb_fpwrap_cdouble_jacobi_theta_1(&cres, cz, ctau, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_jacobi_theta_2(&cres, cz, ctau, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_jacobi_theta_3(&cres, cz, ctau, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_jacobi_theta_4(&cres, cz, ctau, flags));

        CHECK_CDOUBLE(arb_fpwrap_cdouble_dedekind_eta(&cres, ctau, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_modular_j(&cres, ctau, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_modular_lambda(&cres, ctau, flags));
        CHECK_CDOUBLE(arb_fpwrap_cdouble_modular_delta(&cres, ctau, flags));
    }

    TEST_FUNCTION_END(state);
}
