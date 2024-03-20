#include <string.h>
#include <stdlib.h>

/* Include functions *********************************************************/

#include "t-init.c"
#include "t-conversion.c"
#include "t-randtest.c"
#include "t-mul_scalar.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_sparse_vec_init),
    TEST_FUNCTION(gr_sparse_vec_conversion),
    TEST_FUNCTION(gr_sparse_vec_randtest),
    TEST_FUNCTION(gr_sparse_vec_mul_div_scalar),
};

/* main function *************************************************************/

TEST_MAIN(tests)
