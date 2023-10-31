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

#include "t-factor.c"
#include "t-primes_init.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(qsieve_factor),
    TEST_FUNCTION(qsieve_primes_init)
};

/* main function *************************************************************/

TEST_MAIN(tests)
