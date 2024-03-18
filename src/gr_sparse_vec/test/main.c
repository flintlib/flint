#include <string.h>
#include <stdlib.h>

/* Include functions *********************************************************/

#include "t-creation.c"
#include "t-conversion.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_sparse_vec_creation),
    TEST_FUNCTION(gr_sparse_vec_conversion),
};

/* main function *************************************************************/

TEST_MAIN(tests)
