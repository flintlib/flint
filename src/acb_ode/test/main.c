#include "t-apply_diffop.c"
#include "t-sum_divconquer.c"

test_struct tests[] =
{
    TEST_FUNCTION(acb_ode_apply_diffop),
    TEST_FUNCTION(acb_ode_sum_divconquer)
};

TEST_MAIN(tests)
