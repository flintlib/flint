/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>

/* Include functions *********************************************************/

#include "t-convol.c"
#include "t-dft.c"

/* Array of test functions ***************************************************/

int (*test_functions[])(void) =
{
    TEST_FUNCTION(acb_dft_convol),
    TEST_FUNCTION(acb_dft)
};

char acb_dft_convol_name[] = "acb_dft_convol";
char acb_dft_name[] = "acb_dft";

char * test_names[] =
{
    acb_dft_convol_name,
    acb_dft_name
};

/* main function *************************************************************/

#define NUMBER_OF_TESTS (sizeof(test_functions) / sizeof(int (*)(void)))

int
main(int argc, char ** argv)
{
    int ix, jx;

    if (argc < 2)
    {
        for (ix = 0; ix < NUMBER_OF_TESTS; ix++)
            if (test_functions[ix]())
                flint_abort();
    }
    else
    {
        for (ix = 1; ix < argc; ix++)
        {
            for (jx = 0; jx < NUMBER_OF_TESTS; jx++)
            {
                /* If argument equals to test name, run it */
                if (strcmp(argv[ix], test_names[jx]) == 0)
                {
                    if (test_functions[jx]())
                        flint_abort();
                    break;
                }
            }

            if (jx == NUMBER_OF_TESTS)
            {
                fprintf(stderr, "Error: Could not find test function for %s\n", argv[ix]);
                flint_abort();
            }
        }
    }

    return 0;
}
