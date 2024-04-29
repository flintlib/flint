#include <string.h>
#include <stdlib.h>

/* Include functions *********************************************************/

#include "t-init.c"
#include "t-randtest.c"
#include "t-conversion.c"
#include "t-arith.c"
#include "t-mul.c"
#include "t-lu.c"
#include "t-rref.c"
#include "t-solve.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_sparse_mat_init),
    TEST_FUNCTION(gr_sparse_mat_conversion),
    TEST_FUNCTION(gr_sparse_mat_randtest),
    TEST_FUNCTION(gr_sparse_mat_arith),
    TEST_FUNCTION(gr_sparse_mat_mul),
    TEST_FUNCTION(gr_sparse_mat_rref),
    TEST_FUNCTION(gr_sparse_mat_lu),
    TEST_FUNCTION(gr_sparse_mat_solve),
};

/* main function *************************************************************/

TEST_MAIN(tests)
