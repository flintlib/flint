#include <string.h>
#include <stdlib.h>

/* Include functions *********************************************************/

#include "t-init.c"
#include "t-conversion.c"
#include "t-randtest.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_sparse_vec_init),
    TEST_FUNCTION(gr_sparse_vec_conversion),
    TEST_FUNCTION(gr_sparse_vec_randtest),
};

/* main function *************************************************************/

TEST_MAIN(tests)
