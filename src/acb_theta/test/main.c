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

#include "t-agm_hadamard.c"
#include "t-agm_mul.c"
#include "t-agm_mul_tight.c"
#include "t-agm_sqrt.c"
#include "t-all.c"
#include "t-char_dot.c"
#include "t-char_get_a.c"
#include "t-char_is_even.c"
#include "t-char_is_goepel.c"
#include "t-char_is_syzygous.c"
#include "t-dist_a0.c"
#include "t-dist_lat.c"
#include "t-dist_pt.c"
#include "t-eld_border.c"
#include "t-eld_points.c"
#include "t-g2_character.c"
#include "t-g2_chi10.c"
#include "t-g2_chi12.c"
#include "t-g2_chi35.c"
#include "t-g2_chi3_6.c"
#include "t-g2_chi5.c"
#include "t-g2_covariants.c"
#include "t-g2_covariants_lead.c"
#include "t-g2_detk_symj.c"
#include "t-g2_jet_naive_1.c"
#include "t-g2_psi4.c"
#include "t-g2_psi6.c"
#include "t-g2_sextic.c"
#include "t-g2_sextic_chi5.c"
#include "t-g2_transvectant.c"
#include "t-g2_transvectant_lead.c"
#include "t-jet_all.c"
#include "t-jet_compose.c"
#include "t-jet_error_bounds.c"
#include "t-jet_mul.c"
#include "t-jet_naive_00.c"
#include "t-jet_naive_all.c"
#include "t-jet_naive_fixed_ab.c"
#include "t-jet_naive_radius.c"
#include "t-jet_ql_all.c"
#include "t-jet_ql_bounds.c"
#include "t-jet_ql_finite_diff.c"
#include "t-jet_ql_radius.c"
#include "t-jet_tuples.c"
#include "t-naive_00.c"
#include "t-naive_all.c"
#include "t-naive_fixed_ab.c"
#include "t-naive_fixed_a.c"
#include "t-naive_radius.c"
#include "t-naive_reduce.c"
#include "t-naive_term.c"
#include "t-ql_a0.c"
#include "t-ql_a0_split.c"
#include "t-ql_a0_steps.c"
#include "t-ql_all.c"
#include "t-ql_reduce.c"
#include "t-siegel_cocycle.c"
#include "t-siegel_is_reduced.c"
#include "t-siegel_reduce.c"
#include "t-siegel_transform.c"
#include "t-siegel_transform_z.c"
#include "t-sp2gz_decompose.c"
#include "t-sp2gz_inv.c"
#include "t-sp2gz_is_correct.c"
#include "t-sp2gz_set_blocks.c"
#include "t-transform_char.c"
#include "t-transform_kappa.c"
#include "t-transform_proj.c"
#include "t-transform_sqrtdet.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(acb_theta_agm_hadamard),
    TEST_FUNCTION(acb_theta_agm_mul),
    TEST_FUNCTION(acb_theta_agm_mul_tight),
    TEST_FUNCTION(acb_theta_agm_sqrt),
    TEST_FUNCTION(acb_theta_all),
    TEST_FUNCTION(acb_theta_char_dot),
    TEST_FUNCTION(acb_theta_char_get_a),
    TEST_FUNCTION(acb_theta_char_is_even),
    TEST_FUNCTION(acb_theta_char_is_goepel),
    TEST_FUNCTION(acb_theta_char_is_syzygous),
    TEST_FUNCTION(acb_theta_dist_a0),
    TEST_FUNCTION(acb_theta_dist_lat),
    TEST_FUNCTION(acb_theta_dist_pt),
    TEST_FUNCTION(acb_theta_eld_border),
    TEST_FUNCTION(acb_theta_eld_points),
    TEST_FUNCTION(acb_theta_g2_character),
    TEST_FUNCTION(acb_theta_g2_chi10),
    TEST_FUNCTION(acb_theta_g2_chi12),
    TEST_FUNCTION(acb_theta_g2_chi35),
    TEST_FUNCTION(acb_theta_g2_chi3_6),
    TEST_FUNCTION(acb_theta_g2_chi5),
    TEST_FUNCTION(acb_theta_g2_covariants),
    TEST_FUNCTION(acb_theta_g2_covariants_lead),
    TEST_FUNCTION(acb_theta_g2_detk_symj),
    TEST_FUNCTION(acb_theta_g2_jet_naive_1),
    TEST_FUNCTION(acb_theta_g2_psi4),
    TEST_FUNCTION(acb_theta_g2_psi6),
    TEST_FUNCTION(acb_theta_g2_sextic),
    TEST_FUNCTION(acb_theta_g2_sextic_chi5),
    TEST_FUNCTION(acb_theta_g2_transvectant),
    TEST_FUNCTION(acb_theta_g2_transvectant_lead),
    TEST_FUNCTION(acb_theta_jet_all),
    TEST_FUNCTION(acb_theta_jet_compose),
    TEST_FUNCTION(acb_theta_jet_error_bounds),
    TEST_FUNCTION(acb_theta_jet_mul),
    TEST_FUNCTION(acb_theta_jet_naive_00),
    TEST_FUNCTION(acb_theta_jet_naive_all),
    TEST_FUNCTION(acb_theta_jet_naive_fixed_ab),
    TEST_FUNCTION(acb_theta_jet_naive_radius),
    TEST_FUNCTION(acb_theta_jet_ql_all),
    TEST_FUNCTION(acb_theta_jet_ql_bounds),
    TEST_FUNCTION(acb_theta_jet_ql_finite_diff),
    TEST_FUNCTION(acb_theta_jet_ql_radius),
    TEST_FUNCTION(acb_theta_jet_tuples),
    TEST_FUNCTION(acb_theta_naive_00),
    TEST_FUNCTION(acb_theta_naive_all),
    TEST_FUNCTION(acb_theta_naive_fixed_ab),
    TEST_FUNCTION(acb_theta_naive_fixed_a),
    TEST_FUNCTION(acb_theta_naive_radius),
    TEST_FUNCTION(acb_theta_naive_reduce),
    TEST_FUNCTION(acb_theta_naive_term),
    TEST_FUNCTION(acb_theta_ql_a0),
    TEST_FUNCTION(acb_theta_ql_a0_split),
    TEST_FUNCTION(acb_theta_ql_a0_steps),
    TEST_FUNCTION(acb_theta_ql_all),
    TEST_FUNCTION(acb_theta_ql_reduce),
    TEST_FUNCTION(acb_theta_siegel_cocycle),
    TEST_FUNCTION(acb_theta_siegel_is_reduced),
    TEST_FUNCTION(acb_theta_siegel_reduce),
    TEST_FUNCTION(acb_theta_siegel_transform),
    TEST_FUNCTION(acb_theta_siegel_transform_z),
    TEST_FUNCTION(acb_theta_sp2gz_decompose),
    TEST_FUNCTION(acb_theta_sp2gz_inv),
    TEST_FUNCTION(acb_theta_sp2gz_is_correct),
    TEST_FUNCTION(acb_theta_sp2gz_set_blocks),
    TEST_FUNCTION(acb_theta_transform_char),
    TEST_FUNCTION(acb_theta_transform_kappa),
    TEST_FUNCTION(acb_theta_transform_proj),
    TEST_FUNCTION(acb_theta_transform_sqrtdet)
};

/* main function *************************************************************/

TEST_MAIN(tests)
