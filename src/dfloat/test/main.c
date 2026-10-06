/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "dfloat.h"
#include "test_helpers.h"

/* The dfloat arithmetic needs doubles without excess precision; where
   the module was compiled without that (x87), dfloat_is_supported()
   returns 0 and the tests are skipped. */
#define DFLOAT_TEST_SKIP_IF_UNSUPPORTED(state) \
    do { \
        if (!dfloat_is_supported()) \
        { \
            FLINT_TEST_CLEAR(state); \
            printf("%.*s(" _YELLOW_B "SKIPPED" _RESET ")\n", 54, _test_io_string_); \
            return 0; \
        } \
    } while (0)


#include "t-ball_arith.c"
#include "t-canonical.c"
#include "t-complex.c"
#include "t-elementary.c"
#include "t-exp.c"
#include "t-gr.c"
#include "t-gr_vec.c"
#include "t-logatan.c"
#include "t-nonfinite.c"
#include "t-platt_smk.c"
#include "t-powsum.c"
#include "t-predicates.c"
#include "t-special.c"
#include "t-trig.c"
#include "t-vec.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(ball_arith),
    TEST_FUNCTION(canonical),
    TEST_FUNCTION(complex),
    TEST_FUNCTION(elementary),
    TEST_FUNCTION(exp),
    TEST_FUNCTION(gr),
    TEST_FUNCTION(gr_vec),
    TEST_FUNCTION(logatan),
    TEST_FUNCTION(nonfinite),
    TEST_FUNCTION(platt_smk),
    TEST_FUNCTION(powsum),
    TEST_FUNCTION(predicates),
    TEST_FUNCTION(special),
    TEST_FUNCTION(trig),
    TEST_FUNCTION(vec),
};

/* main function *************************************************************/

TEST_MAIN(tests)
