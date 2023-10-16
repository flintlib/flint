/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>

/* Include functions *********************************************************/

#include "t-acos.c"
#include "t-acosh.c"
#include "t-agm1.c"
#include "t-agm.c"
#include "t-approx_dot.c"
#include "t-asin.c"
#include "t-asinh.c"
#include "t-atan.c"
#include "t-atanh.c"
#include "t-barnes_g.c"
#include "t-bernoulli_poly_ui.c"
#include "t-chebyshev_t_ui.c"
#include "t-chebyshev_u_ui.c"
#include "t-cos_pi.c"
#include "t-cot.c"
#include "t-coth.c"
#include "t-cot_pi.c"
#include "t-csc.c"
#include "t-csch.c"
#include "t-csc_pi.c"
#include "t-csgn.c"
#include "t-digamma.c"
#include "t-div.c"
#include "t-dot.c"
#include "t-dot_fmpz.c"
#include "t-dot_si.c"
#include "t-dot_siui.c"
#include "t-dot_ui.c"
#include "t-dot_uiui.c"
#include "t-exp.c"
#include "t-exp_invexp.c"
#include "t-expm1.c"
#include "t-exp_pi_i.c"
#include "t-gamma.c"
#include "t-get_abs_lbound_arf.c"
#include "t-get_abs_ubound_arf.c"
#include "t-get_mag.c"
#include "t-get_mag_lower.c"
#include "t-inv.c"
#include "t-lambertw.c"
#include "t-lgamma.c"
#include "t-log1p.c"
#include "t-log.c"
#include "t-log_sin_pi.c"
#include "t-mul.c"
#include "t-mul_naive.c"
#include "t-polygamma.c"
#include "t-pow.c"
#include "t-pow_fmpz.c"
#include "t-quadratic_roots_fmpz.c"
#include "t-rel_accuracy_bits.c"
#include "t-rgamma.c"
#include "t-rising2_ui.c"
#include "t-rising_ui.c"
#include "t-rising_ui_get_mag.c"
#include "t-root_ui.c"
#include "t-rsqrt.c"
#include "t-sec.c"
#include "t-sech.c"
#include "t-sgn.c"
#include "t-sinc.c"
#include "t-sin_cos.c"
#include "t-sinc_pi.c"
#include "t-sinh_cosh.c"
#include "t-sin_pi.c"
#include "t-sqrt.c"
#include "t-tan.c"
#include "t-tanh.c"
#include "t-tan_pi.c"
#include "t-vec_unit_roots.c"
#include "t-zeta.c"

/* Array of test functions ***************************************************/

int (*test_functions[])(void) =
{
    TEST_FUNCTION(acb_acos),
    TEST_FUNCTION(acb_acosh),
    TEST_FUNCTION(acb_agm1),
    TEST_FUNCTION(acb_agm),
    TEST_FUNCTION(acb_approx_dot),
    TEST_FUNCTION(acb_asin),
    TEST_FUNCTION(acb_asinh),
    TEST_FUNCTION(acb_atan),
    TEST_FUNCTION(acb_atanh),
    TEST_FUNCTION(acb_barnes_g),
    TEST_FUNCTION(acb_bernoulli_poly_ui),
    TEST_FUNCTION(acb_chebyshev_t_ui),
    TEST_FUNCTION(acb_chebyshev_u_ui),
    TEST_FUNCTION(acb_cos_pi),
    TEST_FUNCTION(acb_cot),
    TEST_FUNCTION(acb_coth),
    TEST_FUNCTION(acb_cot_pi),
    TEST_FUNCTION(acb_csc),
    TEST_FUNCTION(acb_csch),
    TEST_FUNCTION(acb_csc_pi),
    TEST_FUNCTION(acb_csgn),
    TEST_FUNCTION(acb_digamma),
    TEST_FUNCTION(acb_div),
    TEST_FUNCTION(acb_dot),
    TEST_FUNCTION(acb_dot_fmpz),
    TEST_FUNCTION(acb_dot_si),
    TEST_FUNCTION(acb_dot_siui),
    TEST_FUNCTION(acb_dot_ui),
    TEST_FUNCTION(acb_dot_uiui),
    TEST_FUNCTION(acb_exp),
    TEST_FUNCTION(acb_exp_invexp),
    TEST_FUNCTION(acb_expm1),
    TEST_FUNCTION(acb_exp_pi_i),
    TEST_FUNCTION(acb_gamma),
    TEST_FUNCTION(acb_get_abs_lbound_arf),
    TEST_FUNCTION(acb_get_abs_ubound_arf),
    TEST_FUNCTION(acb_get_mag),
    TEST_FUNCTION(acb_get_mag_lower),
    TEST_FUNCTION(acb_inv),
    TEST_FUNCTION(acb_lambertw),
    TEST_FUNCTION(acb_lgamma),
    TEST_FUNCTION(acb_log1p),
    TEST_FUNCTION(acb_log),
    TEST_FUNCTION(acb_log_sin_pi),
    TEST_FUNCTION(acb_mul),
    TEST_FUNCTION(acb_mul_naive),
    TEST_FUNCTION(acb_polygamma),
    TEST_FUNCTION(acb_pow),
    TEST_FUNCTION(acb_pow_fmpz),
    TEST_FUNCTION(acb_quadratic_roots_fmpz),
    TEST_FUNCTION(acb_rel_accuracy_bits),
    TEST_FUNCTION(acb_rgamma),
    TEST_FUNCTION(acb_rising2_ui),
    TEST_FUNCTION(acb_rising_ui),
    TEST_FUNCTION(acb_rising_ui_get_mag),
    TEST_FUNCTION(acb_root_ui),
    TEST_FUNCTION(acb_rsqrt),
    TEST_FUNCTION(acb_sec),
    TEST_FUNCTION(acb_sech),
    TEST_FUNCTION(acb_sgn),
    TEST_FUNCTION(acb_sinc),
    TEST_FUNCTION(acb_sin_cos),
    TEST_FUNCTION(acb_sinc_pi),
    TEST_FUNCTION(acb_sinh_cosh),
    TEST_FUNCTION(acb_sin_pi),
    TEST_FUNCTION(acb_sqrt),
    TEST_FUNCTION(acb_tan),
    TEST_FUNCTION(acb_tanh),
    TEST_FUNCTION(acb_tan_pi),
    TEST_FUNCTION(acb_vec_unit_roots),
    TEST_FUNCTION(acb_zeta)
};

char acb_acos_name[] = "acb_acos";
char acb_acosh_name[] = "acb_acosh";
char acb_agm1_name[] = "acb_agm1";
char acb_agm_name[] = "acb_agm";
char acb_approx_dot_name[] = "acb_approx_dot";
char acb_asin_name[] = "acb_asin";
char acb_asinh_name[] = "acb_asinh";
char acb_atan_name[] = "acb_atan";
char acb_atanh_name[] = "acb_atanh";
char acb_barnes_g_name[] = "acb_barnes_g";
char acb_bernoulli_poly_ui_name[] = "acb_bernoulli_poly_ui";
char acb_chebyshev_t_ui_name[] = "acb_chebyshev_t_ui";
char acb_chebyshev_u_ui_name[] = "acb_chebyshev_u_ui";
char acb_cos_pi_name[] = "acb_cos_pi";
char acb_cot_name[] = "acb_cot";
char acb_coth_name[] = "acb_coth";
char acb_cot_pi_name[] = "acb_cot_pi";
char acb_csc_name[] = "acb_csc";
char acb_csch_name[] = "acb_csch";
char acb_csc_pi_name[] = "acb_csc_pi";
char acb_csgn_name[] = "acb_csgn";
char acb_digamma_name[] = "acb_digamma";
char acb_div_name[] = "acb_div";
char acb_dot_name[] = "acb_dot";
char acb_dot_fmpz_name[] = "acb_dot_fmpz";
char acb_dot_si_name[] = "acb_dot_si";
char acb_dot_siui_name[] = "acb_dot_siui";
char acb_dot_ui_name[] = "acb_dot_ui";
char acb_dot_uiui_name[] = "acb_dot_uiui";
char acb_exp_name[] = "acb_exp";
char acb_exp_invexp_name[] = "acb_exp_invexp";
char acb_expm1_name[] = "acb_expm1";
char acb_exp_pi_i_name[] = "acb_exp_pi_i";
char acb_gamma_name[] = "acb_gamma";
char acb_get_abs_lbound_arf_name[] = "acb_get_abs_lbound_arf";
char acb_get_abs_ubound_arf_name[] = "acb_get_abs_ubound_arf";
char acb_get_mag_name[] = "acb_get_mag";
char acb_get_mag_lower_name[] = "acb_get_mag_lower";
char acb_inv_name[] = "acb_inv";
char acb_lambertw_name[] = "acb_lambertw";
char acb_lgamma_name[] = "acb_lgamma";
char acb_log1p_name[] = "acb_log1p";
char acb_log_name[] = "acb_log";
char acb_log_sin_pi_name[] = "acb_log_sin_pi";
char acb_mul_name[] = "acb_mul";
char acb_mul_naive_name[] = "acb_mul_naive";
char acb_polygamma_name[] = "acb_polygamma";
char acb_pow_name[] = "acb_pow";
char acb_pow_fmpz_name[] = "acb_pow_fmpz";
char acb_quadratic_roots_fmpz_name[] = "acb_quadratic_roots_fmpz";
char acb_rel_accuracy_bits_name[] = "acb_rel_accuracy_bits";
char acb_rgamma_name[] = "acb_rgamma";
char acb_rising2_ui_name[] = "acb_rising2_ui";
char acb_rising_ui_name[] = "acb_rising_ui";
char acb_rising_ui_get_mag_name[] = "acb_rising_ui_get_mag";
char acb_root_ui_name[] = "acb_root_ui";
char acb_rsqrt_name[] = "acb_rsqrt";
char acb_sec_name[] = "acb_sec";
char acb_sech_name[] = "acb_sech";
char acb_sgn_name[] = "acb_sgn";
char acb_sinc_name[] = "acb_sinc";
char acb_sin_cos_name[] = "acb_sin_cos";
char acb_sinc_pi_name[] = "acb_sinc_pi";
char acb_sinh_cosh_name[] = "acb_sinh_cosh";
char acb_sin_pi_name[] = "acb_sin_pi";
char acb_sqrt_name[] = "acb_sqrt";
char acb_tan_name[] = "acb_tan";
char acb_tanh_name[] = "acb_tanh";
char acb_tan_pi_name[] = "acb_tan_pi";
char acb_vec_unit_roots_name[] = "acb_vec_unit_roots";
char acb_zeta_name[] = "acb_zeta";

char * test_names[] =
{
    acb_acos_name,
    acb_acosh_name,
    acb_agm1_name,
    acb_agm_name,
    acb_approx_dot_name,
    acb_asin_name,
    acb_asinh_name,
    acb_atan_name,
    acb_atanh_name,
    acb_barnes_g_name,
    acb_bernoulli_poly_ui_name,
    acb_chebyshev_t_ui_name,
    acb_chebyshev_u_ui_name,
    acb_cos_pi_name,
    acb_cot_name,
    acb_coth_name,
    acb_cot_pi_name,
    acb_csc_name,
    acb_csch_name,
    acb_csc_pi_name,
    acb_csgn_name,
    acb_digamma_name,
    acb_div_name,
    acb_dot_name,
    acb_dot_fmpz_name,
    acb_dot_si_name,
    acb_dot_siui_name,
    acb_dot_ui_name,
    acb_dot_uiui_name,
    acb_exp_name,
    acb_exp_invexp_name,
    acb_expm1_name,
    acb_exp_pi_i_name,
    acb_gamma_name,
    acb_get_abs_lbound_arf_name,
    acb_get_abs_ubound_arf_name,
    acb_get_mag_name,
    acb_get_mag_lower_name,
    acb_inv_name,
    acb_lambertw_name,
    acb_lgamma_name,
    acb_log1p_name,
    acb_log_name,
    acb_log_sin_pi_name,
    acb_mul_name,
    acb_mul_naive_name,
    acb_polygamma_name,
    acb_pow_name,
    acb_pow_fmpz_name,
    acb_quadratic_roots_fmpz_name,
    acb_rel_accuracy_bits_name,
    acb_rgamma_name,
    acb_rising2_ui_name,
    acb_rising_ui_name,
    acb_rising_ui_get_mag_name,
    acb_root_ui_name,
    acb_rsqrt_name,
    acb_sec_name,
    acb_sech_name,
    acb_sgn_name,
    acb_sinc_name,
    acb_sin_cos_name,
    acb_sinc_pi_name,
    acb_sinh_cosh_name,
    acb_sin_pi_name,
    acb_sqrt_name,
    acb_tan_name,
    acb_tanh_name,
    acb_tan_pi_name,
    acb_vec_unit_roots_name,
    acb_zeta_name
};

/* main function *************************************************************/

#define NUMBER_OF_TESTS (sizeof(test_functions) / sizeof(int (*)(void)))

int
main(int argc, char ** argv)
{
    int ix, jx;

    if (argc < 2)
    {
        for (ix = 0; ix < NUMBER_OF_TESTS; ix++)
            if (test_functions[ix]())
                flint_abort();
    }
    else
    {
        for (ix = 1; ix < argc; ix++)
        {
            for (jx = 0; jx < NUMBER_OF_TESTS; jx++)
            {
                /* If argument equals to test name, run it */
                if (strcmp(argv[ix], test_names[jx]) == 0)
                {
                    if (test_functions[jx]())
                        flint_abort();
                    break;
                }
            }

            if (jx == NUMBER_OF_TESTS)
            {
                fprintf(stderr, "Error: Could not find test function for %s\n", argv[ix]);
                flint_abort();
            }
        }
    }

    return 0;
}
