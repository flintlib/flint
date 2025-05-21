/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-agm_mul.c"
#include "t-agm_mul_tight.c"
#include "t-agm_sqrt.c"
#include "t-char_dot.c"
#include "t-char_set_slong_vec.c"
#include "t-char_shuffle.c"
#include "t-char_table.c"
#include "t-ctx_exp_inv.c"
#include "t-ctx_sqr_inv.c"
#include "t-ctx_tau_dupl.c"
#include "t-ctx_z_add_real.c"
#include "t-ctx_z_dupl.c"
#include "t-eld_border.c"
#include "t-eld_distances.c"
#include "t-eld_points.c"
#include "t-g2_character.c"
#include "t-g2_chi3_6.c"
#include "t-g2_chi5.c"
#include "t-g2_chi35.c"
#include "t-g2_covariants.c"
#include "t-g2_detk_symj.c"
#include "t-g2_even_weight.c"
#include "t-g2_sextic_chi5.c"
#include "t-g2_transvectant.c"
#include "t-jet.c"
#include "t-jet_compose.c"
#include "t-jet_exp_pi_i.c"
#include "t-jet_exp_qf.c"
#include "t-jet_mul.c"
#include "t-jet_notransform.c"
#include "t-jet_tuples.c"
#include "t-ql_exact.c"
#include "t-ql_jet_error.c"
#include "t-ql_jet_fd.c"
#include "t-ql_local_bound.c"
#include "t-ql_lower_dim.c"
#include "t-ql_setup.c"
#include "t-reduce_z.c"
#include "t-siegel_cocycle.c"
#include "t-siegel_is_reduced.c"
#include "t-siegel_kappa.c"
#include "t-siegel_reduce.c"
#include "t-siegel_transform.c"
#include "t-sp2gz_decompose.c"
#include "t-sp2gz_inv.c"
#include "t-sp2gz_is_correct.c"
#include "t-sp2gz_set_blocks.c"
#include "t-sum.c"
#include "t-sum_jet.c"
#include "t-sum_jet_radius.c"
#include "t-sum_radius.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(acb_theta_agm_mul),
    TEST_FUNCTION(acb_theta_agm_mul_tight),
    TEST_FUNCTION(acb_theta_agm_sqrt),
    TEST_FUNCTION(acb_theta_char_dot),
    TEST_FUNCTION(acb_theta_char_set_slong_vec),
    TEST_FUNCTION(acb_theta_char_shuffle),
    TEST_FUNCTION(acb_theta_char_table),
    TEST_FUNCTION(acb_theta_ctx_exp_inv),
    TEST_FUNCTION(acb_theta_ctx_sqr_inv),
    TEST_FUNCTION(acb_theta_ctx_tau_dupl),
    TEST_FUNCTION(acb_theta_ctx_z_add_real),
    TEST_FUNCTION(acb_theta_ctx_z_dupl),
    TEST_FUNCTION(acb_theta_eld_border),
    TEST_FUNCTION(acb_theta_eld_distances),
    TEST_FUNCTION(acb_theta_eld_points),
    TEST_FUNCTION(acb_theta_g2_character),
    TEST_FUNCTION(acb_theta_g2_chi3_6),
    TEST_FUNCTION(acb_theta_g2_chi5),
    TEST_FUNCTION(acb_theta_g2_chi35),
    TEST_FUNCTION(acb_theta_g2_covariants),
    TEST_FUNCTION(acb_theta_g2_detk_symj),
    TEST_FUNCTION(acb_theta_g2_even_weight),
    TEST_FUNCTION(acb_theta_g2_sextic_chi5),
    TEST_FUNCTION(acb_theta_g2_transvectant),
    TEST_FUNCTION(acb_theta_jet),
    TEST_FUNCTION(acb_theta_jet_compose),
    TEST_FUNCTION(acb_theta_jet_exp_pi_i),
    TEST_FUNCTION(acb_theta_jet_exp_qf),
    TEST_FUNCTION(acb_theta_jet_mul),
    TEST_FUNCTION(acb_theta_jet_notransform),
    TEST_FUNCTION(acb_theta_jet_tuples),
    TEST_FUNCTION(acb_theta_ql_exact),
    TEST_FUNCTION(acb_theta_ql_jet_error),
    TEST_FUNCTION(acb_theta_ql_jet_fd),
    TEST_FUNCTION(acb_theta_ql_local_bound),
    TEST_FUNCTION(acb_theta_ql_lower_dim),
    TEST_FUNCTION(acb_theta_ql_setup),
    TEST_FUNCTION(acb_theta_reduce_z),
    TEST_FUNCTION(acb_theta_siegel_cocycle),
    TEST_FUNCTION(acb_theta_siegel_is_reduced),
    TEST_FUNCTION(acb_theta_siegel_kappa),
    TEST_FUNCTION(acb_theta_siegel_reduce),
    TEST_FUNCTION(acb_theta_siegel_transform),
    TEST_FUNCTION(acb_theta_sp2gz_decompose),
    TEST_FUNCTION(acb_theta_sp2gz_inv),
    TEST_FUNCTION(acb_theta_sp2gz_is_correct),
    TEST_FUNCTION(acb_theta_sp2gz_set_blocks),
    TEST_FUNCTION(acb_theta_sum),
    TEST_FUNCTION(acb_theta_sum_jet),
    TEST_FUNCTION(acb_theta_sum_jet_radius),
    TEST_FUNCTION(acb_theta_sum_radius),
};

/* main function *************************************************************/

TEST_MAIN(tests)
