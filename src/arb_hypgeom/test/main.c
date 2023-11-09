/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdlib.h>

/* Include functions *********************************************************/

#include "t-1f1_integration.c"
#include "t-2f1_integration.c"
#include "t-airy_zero.c"
#include "t-bessel_i_integration.c"
#include "t-bessel_j.c"
#include "t-bessel_k_integration.c"
#include "t-bessel_y.c"
#include "t-central_bin_ui.c"
#include "t-ci.c"
#include "t-coulomb_series.c"
#include "t-erf.c"
#include "t-erfinv.c"
#include "t-gamma_fmpq.c"
#include "t-gamma_lower_sum_rs.c"
#include "t-gamma_stirling_sum.c"
#include "t-gamma_taylor.c"
#include "t-gamma_taylor_tab.c"
#include "t-gamma_upper_fmpq.c"
#include "t-gamma_upper_integration.c"
#include "t-gamma_upper_sum_rs.c"
#include "t-legendre_p_ui_asymp.c"
#include "t-legendre_p_ui.c"
#include "t-legendre_p_ui_deriv_bound.c"
#include "t-legendre_p_ui_one.c"
#include "t-legendre_p_ui_rec.c"
#include "t-legendre_p_ui_root.c"
#include "t-legendre_p_ui_zero.c"
#include "t-lgamma.c"
#include "t-rising_ui.c"
#include "t-rising_ui_jet.c"
#include "t-si.c"
#include "t-sum_fmpq_arb.c"
#include "t-sum_fmpq_imag_arb.c"
#include "t-u_integration.c"
#include "t-wrappers.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(arb_hypgeom_1f1_integration),
    TEST_FUNCTION(arb_hypgeom_2f1_integration),
    TEST_FUNCTION(arb_hypgeom_airy_zero),
    TEST_FUNCTION(arb_hypgeom_bessel_i_integration),
    TEST_FUNCTION(arb_hypgeom_bessel_j),
    TEST_FUNCTION(arb_hypgeom_bessel_k_integration),
    TEST_FUNCTION(arb_hypgeom_bessel_y),
    TEST_FUNCTION(arb_hypgeom_central_bin_ui),
    TEST_FUNCTION(arb_hypgeom_ci),
    TEST_FUNCTION(arb_hypgeom_coulomb_series),
    TEST_FUNCTION(arb_hypgeom_erf),
    TEST_FUNCTION(arb_hypgeom_erfinv),
    TEST_FUNCTION(arb_hypgeom_gamma_fmpq),
    TEST_FUNCTION(arb_hypgeom_gamma_lower_sum_rs),
    TEST_FUNCTION(arb_hypgeom_gamma_stirling_sum),
    TEST_FUNCTION(arb_hypgeom_gamma_taylor),
    TEST_FUNCTION(arb_hypgeom_gamma_taylor_tab),
    TEST_FUNCTION(arb_hypgeom_gamma_upper_fmpq),
    TEST_FUNCTION(arb_hypgeom_gamma_upper_integration),
    TEST_FUNCTION(arb_hypgeom_gamma_upper_sum_rs),
    TEST_FUNCTION(arb_hypgeom_legendre_p_ui_asymp),
    TEST_FUNCTION(arb_hypgeom_legendre_p_ui),
    TEST_FUNCTION(arb_hypgeom_legendre_p_ui_deriv_bound),
    TEST_FUNCTION(arb_hypgeom_legendre_p_ui_one),
    TEST_FUNCTION(arb_hypgeom_legendre_p_ui_rec),
    TEST_FUNCTION(arb_hypgeom_legendre_p_ui_root),
    TEST_FUNCTION(arb_hypgeom_legendre_p_ui_zero),
    TEST_FUNCTION(arb_hypgeom_lgamma),
    TEST_FUNCTION(arb_hypgeom_rising_ui),
    TEST_FUNCTION(arb_hypgeom_rising_ui_jet),
    TEST_FUNCTION(arb_hypgeom_si),
    TEST_FUNCTION(arb_hypgeom_sum_fmpq_arb),
    TEST_FUNCTION(arb_hypgeom_sum_fmpq_imag_arb),
    TEST_FUNCTION(arb_hypgeom_u_integration),
    TEST_FUNCTION(arb_hypgeom_wrappers)
};

/* main function *************************************************************/

TEST_MAIN(tests)
