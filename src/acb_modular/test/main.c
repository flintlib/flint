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

#include "t-delta.c"
#include "t-eisenstein.c"
#include "t-elliptic_e.c"
#include "t-elliptic_k.c"
#include "t-elliptic_p.c"
#include "t-elliptic_p_zpx.c"
#include "t-epsilon_arg.c"
#include "t-eta.c"
#include "t-fundamental_domain_approx.c"
#include "t-hilbert_class_poly.c"
#include "t-j.c"
#include "t-lambda.c"
#include "t-psl2z_inv.c"
#include "t-psl2z_mul.c"
#include "t-theta.c"
#include "t-theta_const_sum_rs.c"
#include "t-theta_jet.c"
#include "t-theta_series.c"
#include "t-theta_sum.c"
#include "t-transform.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(acb_modular_delta),
    TEST_FUNCTION(acb_modular_eisenstein),
    TEST_FUNCTION(acb_modular_elliptic_e),
    TEST_FUNCTION(acb_modular_elliptic_k),
    TEST_FUNCTION(acb_modular_elliptic_p),
    TEST_FUNCTION(acb_modular_elliptic_p_zpx),
    TEST_FUNCTION(acb_modular_epsilon_arg),
    TEST_FUNCTION(acb_modular_eta),
    TEST_FUNCTION(acb_modular_fundamental_domain_approx),
    TEST_FUNCTION(acb_modular_hilbert_class_poly),
    TEST_FUNCTION(acb_modular_j),
    TEST_FUNCTION(acb_modular_lambda),
    TEST_FUNCTION(acb_modular_psl2z_inv),
    TEST_FUNCTION(acb_modular_psl2z_mul),
    TEST_FUNCTION(acb_modular_theta),
    TEST_FUNCTION(acb_modular_theta_const_sum_rs),
    TEST_FUNCTION(acb_modular_theta_jet),
    TEST_FUNCTION(acb_modular_theta_series),
    TEST_FUNCTION(acb_modular_theta_sum),
    TEST_FUNCTION(acb_modular_transform)
};

/* main function *************************************************************/

TEST_MAIN(tests)
