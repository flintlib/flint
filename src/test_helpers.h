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

#define TEST_FUNCTION(label) CAT(test, label)

#define TEST_FUNCTION_START(label)                      \
int CAT(test, label)(void)                              \
{                                                       \
    FLINT_TEST_INIT(state);                             \
    printf(#label "....");                              \
    fflush(stdout);                                     \

#define TEST_FUNCTION_END                               \
    FLINT_TEST_CLEANUP(state);                          \
    printf("PASS\n");                                   \
    return 0;                                           \
}

#define TEST_FUNCTION_END_SKIPPED                       \
    FLINT_TEST_CLEANUP(state);                          \
    printf("SKIPPED\n");                                \
    return 0;                                           \
}

#endif
