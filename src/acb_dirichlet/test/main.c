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

#include "t-backlund_s_bound.c"
#include "t-backlund_s.c"
#include "t-backlund_s_gram.c"
#include "t-chi.c"
#include "t-dft.c"
#include "t-eta.c"
#include "t-euler_product_real_ui.c"
#include "t-gauss.c"
#include "t-gram_point.c"
#include "t-hardy_theta_series.c"
#include "t-hardy_z.c"
#include "t-hardy_z_series.c"
#include "t-hardy_z_zero.c"
#include "t-hardy_z_zeros.c"
#include "t-hurwitz.c"
#include "t-hurwitz_precomp.c"
#include "t-isolate_hardy_z_zero.c"
#include "t-jacobi.c"
#include "t-l.c"
#include "t-lerch_phi.c"
#include "t-l_euler_product.c"
#include "t-l_fmpq_afe.c"
#include "t-l_fmpq.c"
#include "t-l_hurwitz.c"
#include "t-l_jet.c"
#include "t-l_series.c"
#include "t-l_vec_hurwitz.c"
#include "t-platt_beta.c"
#include "t-platt_hardy_z_zeros.c"
#include "t-platt_local_hardy_z_zeros.c"
#include "t-platt_multieval.c"
#include "t-platt_multieval_threaded.c"
#include "t-platt_ws_interpolation.c"
#include "t-platt_zeta_zeros.c"
#include "t-powsum_smooth.c"
#include "t-roots.c"
#include "t-stieltjes.c"
#include "t-thetanull.c"
#include "t-turing_method_bound.c"
#include "t-xi.c"
#include "t-zeta_bound.c"
#include "t-zeta_jet_rs.c"
#include "t-zeta_nzeros.c"
#include "t-zeta_nzeros_gram.c"
#include "t-zeta_rs.c"
#include "t-zeta_rs_r.c"
#include "t-zeta_zero.c"
#include "t-zeta_zeros.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(acb_dirichlet_backlund_s_bound),
    TEST_FUNCTION(acb_dirichlet_backlund_s),
    TEST_FUNCTION(acb_dirichlet_backlund_s_gram),
    TEST_FUNCTION(acb_dirichlet_chi),
    TEST_FUNCTION(acb_dirichlet_dft),
    TEST_FUNCTION(acb_dirichlet_eta),
    TEST_FUNCTION(acb_dirichlet_euler_product_real_ui),
    TEST_FUNCTION(acb_dirichlet_gauss),
    TEST_FUNCTION(acb_dirichlet_gram_point),
    TEST_FUNCTION(acb_dirichlet_hardy_theta_series),
    TEST_FUNCTION(acb_dirichlet_hardy_z),
    TEST_FUNCTION(acb_dirichlet_hardy_z_series),
    TEST_FUNCTION(acb_dirichlet_hardy_z_zero),
    TEST_FUNCTION(acb_dirichlet_hardy_z_zeros),
    TEST_FUNCTION(acb_dirichlet_hurwitz),
    TEST_FUNCTION(acb_dirichlet_hurwitz_precomp),
    TEST_FUNCTION(acb_dirichlet_isolate_hardy_z_zero),
    TEST_FUNCTION(acb_dirichlet_jacobi),
    TEST_FUNCTION(acb_dirichlet_l),
    TEST_FUNCTION(acb_dirichlet_lerch_phi),
    TEST_FUNCTION(acb_dirichlet_l_euler_product),
    TEST_FUNCTION(acb_dirichlet_l_fmpq_afe),
    TEST_FUNCTION(acb_dirichlet_l_fmpq),
    TEST_FUNCTION(acb_dirichlet_l_hurwitz),
    TEST_FUNCTION(acb_dirichlet_l_jet),
    TEST_FUNCTION(acb_dirichlet_l_series),
    TEST_FUNCTION(acb_dirichlet_l_vec_hurwitz),
    TEST_FUNCTION(acb_dirichlet_platt_beta),
    TEST_FUNCTION(acb_dirichlet_platt_hardy_z_zeros),
    TEST_FUNCTION(acb_dirichlet_platt_local_hardy_z_zeros),
    TEST_FUNCTION(acb_dirichlet_platt_multieval),
    TEST_FUNCTION(acb_dirichlet_platt_multieval_threaded),
    TEST_FUNCTION(acb_dirichlet_platt_ws_interpolation),
    TEST_FUNCTION(acb_dirichlet_platt_zeta_zeros),
    TEST_FUNCTION(acb_dirichlet_powsum_smooth),
    TEST_FUNCTION(acb_dirichlet_roots),
    TEST_FUNCTION(acb_dirichlet_stieltjes),
    TEST_FUNCTION(acb_dirichlet_thetanull),
    TEST_FUNCTION(acb_dirichlet_turing_method_bound),
    TEST_FUNCTION(acb_dirichlet_xi),
    TEST_FUNCTION(acb_dirichlet_zeta_bound),
    TEST_FUNCTION(acb_dirichlet_zeta_jet_rs),
    TEST_FUNCTION(acb_dirichlet_zeta_nzeros),
    TEST_FUNCTION(acb_dirichlet_zeta_nzeros_gram),
    TEST_FUNCTION(acb_dirichlet_zeta_rs),
    TEST_FUNCTION(acb_dirichlet_zeta_rs_r),
    TEST_FUNCTION(acb_dirichlet_zeta_zero),
    TEST_FUNCTION(acb_dirichlet_zeta_zeros)
};

/* main function *************************************************************/

TEST_MAIN(tests)
