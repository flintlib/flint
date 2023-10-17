/*
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H

#include "templates.h"

/* easy way to test a condition in test code */
#define FLINT_TEST(e)                                   \
    do {                                                \
        if (!(e))                                       \
        {                                               \
            flint_printf("test %s FAILED\n", #e);       \
            flint_printf("%s:%d\n", __FILE__, __LINE__);\
            flint_abort();                              \
        }                                               \
    } while (0)

/* test function macro *******************************************************/

typedef struct
{
    int (* test_function)(void);
    char * name;
}
test_struct;

#define TEST_FUNCTION(label) { CAT(test, label), TEMPLATE_STR(label) }

#define TEST_FUNCTION_START(label, state)               \
int CAT(test, label)(void)                              \
{                                                       \
    FLINT_TEST_INIT(state);                             \
    printf(#label "....");                              \
    fflush(stdout);                                     \

#define TEST_FUNCTION_END(state)                        \
    FLINT_TEST_CLEANUP(state);                          \
    printf("PASS\n");                                   \
    return 0;                                           \
}

#define TEST_FUNCTION_END_SKIPPED(state)                \
    FLINT_TEST_CLEANUP(state);                          \
    printf("SKIPPED\n");                                \
    return 0;                                           \
}

#define TEST_MAIN(tests)                                                    \
int                                                                         \
main(int argc, char ** argv)                                                \
{                                                                           \
    int ix, jx;                                                             \
                                                                            \
    /* If no arguments where put in, run them all. Else, check which */     \
    /* functions should be runned. */                                       \
    if (argc < 2)                                                           \
    {                                                                       \
        for (ix = 0; ix < (sizeof(tests) / sizeof(test_struct)); ix++)      \
            if ((tests)[ix].test_function())                                \
                flint_abort();                                              \
    }                                                                       \
    else                                                                    \
    {                                                                       \
        for (jx = 1; jx < argc; jx++)                                       \
        {                                                                   \
            for (ix = 0; ix < (sizeof(tests) / sizeof(test_struct)); ix++)  \
            {                                                               \
                /* If argument equals to test name, run it */               \
                if (strcmp(argv[jx], (tests)[ix].name) == 0)                \
                {                                                           \
                    if ((tests)[ix].test_function())                        \
                        flint_abort();                                      \
                    break;                                                  \
                }                                                           \
            }                                                               \
                                                                            \
            if (ix == (sizeof(tests) / sizeof(test_struct)))                \
            {                                                               \
                fprintf(stderr,                                             \
                        "Error: Could not find test function for %s\n",     \
                        argv[jx]);                                          \
                flint_abort();                                              \
            }                                                               \
        }                                                                   \
    }                                                                       \
                                                                            \
    return 0;                                                               \
}

#endif
