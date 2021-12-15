/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H

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

#endif

