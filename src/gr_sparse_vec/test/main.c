#include <string.h>
#include <stdlib.h>

/* Include functions *********************************************************/

#include "t-creation.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_sparse_vec_creation),
};

/* main function *************************************************************/

TEST_MAIN(tests)