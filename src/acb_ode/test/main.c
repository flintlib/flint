#include "t-apply_diffop.c"
#include "t-exponents.c"
#include "t-exponents_ordinary.c"
#include "t-exponents_randtest.c"
#include "t-fundamental_matrix.c"
#include "t-bound_rat.c"
#include "t-indicial_polynomial_from_exponents.c"
#include "t-poly_negdivrevhigh.c"
#include "t-poly_taylor_shift_aps_trunc.c"
#include "t-sum_divconquer.c"

test_struct tests[] =
{
    TEST_FUNCTION(acb_ode_apply_diffop),
    TEST_FUNCTION(acb_ode_bound_rat),
    TEST_FUNCTION(acb_ode_exponents),
    TEST_FUNCTION(acb_ode_exponents_ordinary),
    TEST_FUNCTION(acb_ode_exponents_randtest),
    TEST_FUNCTION(acb_ode_fundamental_matrix),
    TEST_FUNCTION(acb_ode_indicial_polynomial_from_exponents),
    TEST_FUNCTION(acb_ode_poly_negdivrevhigh),
    TEST_FUNCTION(acb_ode_poly_taylor_shift_aps_trunc)
    // TEST_FUNCTION(acb_ode_sum_divconquer)
};

TEST_MAIN(tests)
