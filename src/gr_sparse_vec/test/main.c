#include <string.h>
#include <stdlib.h>

/* Include functions *********************************************************/

#include "t-arith.c"
#include "t-conversion.c"
#include "t-dot.c"
#include "t-init.c"
#include "t-randtest.c"
#include "t-sum-prod.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_sparse_vec_init),
    TEST_FUNCTION(gr_sparse_vec_conversion),
    TEST_FUNCTION(gr_sparse_vec_randtest),
    TEST_FUNCTION(gr_sparse_vec_arith),
    TEST_FUNCTION(gr_sparse_vec_dot),
    TEST_FUNCTION(gr_sparse_vec_sum_prod),
};

/* main function *************************************************************/

TEST_MAIN(tests)
