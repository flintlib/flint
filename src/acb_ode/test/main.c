/*
    Copyright (C) 2025, 2026 Marc Mezzarobba

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "t-apply_diffop.c"
#include "t-bound_rat.c"
#include "t-exponents.c"
#include "t-exponents_ordinary.c"
#include "t-exponents_randtest.c"
#include "t-fundamental_matrix.c"
#include "t-fundamental_matrix_convergence.c"
#include "t-fundamental_matrix_vs_ref.c"
#include "t-ind_lbound_eval.c"
#include "t-indicial_polynomial_from_exponents.c"
#include "t-lc_roots.c"
#include "t-local_diffop.c"
#include "t-poly_negdivrevhigh.c"
#include "t-poly_taylor_shift_aps_trunc.c"
#include "t-randtest_acb.c"
#include "t-solution_growth.c"
#include "t-solution_series.c"
#include "t-sum_divconquer.c"
#include "t-sum_forward.c"
#include "t-sum_update_tail_bounds.c"
#include "t-transition_matrix_fixed.c"

test_struct tests[] =
{
    TEST_FUNCTION(acb_ode_apply_diffop),
    TEST_FUNCTION(acb_ode_bound_rat),
    TEST_FUNCTION(acb_ode_exponents),
    TEST_FUNCTION(acb_ode_exponents_ordinary),
    TEST_FUNCTION(acb_ode_exponents_randtest),
    TEST_FUNCTION(acb_ode_fundamental_matrix),
    TEST_FUNCTION(acb_ode_fundamental_matrix_convergence),
    TEST_FUNCTION(acb_ode_fundamental_matrix_vs_ref),
    TEST_FUNCTION(acb_ode_ind_lbound_eval),
    TEST_FUNCTION(acb_ode_indicial_polynomial_from_exponents),
    TEST_FUNCTION(acb_ode_lc_roots),
    TEST_FUNCTION(acb_ode_local_diffop),
    TEST_FUNCTION(acb_ode_poly_negdivrevhigh),
    TEST_FUNCTION(acb_ode_poly_taylor_shift_aps_trunc),
    TEST_FUNCTION(acb_ode_randtest_acb),
    TEST_FUNCTION(acb_ode_solution_series),
    TEST_FUNCTION(acb_ode_solution_growth),
    TEST_FUNCTION(acb_ode_sum_divconquer),
    TEST_FUNCTION(acb_ode_sum_forward),
    TEST_FUNCTION(acb_ode_sum_update_tail_bounds),
    TEST_FUNCTION(acb_ode_transition_matrix_fixed)
};

TEST_MAIN(tests)
