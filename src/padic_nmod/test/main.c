/*
    Copyright (C) 2026 Rubén Muñoz--Bertrand

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "t-add.c"
#include "t-div.c"
#include "t-inv.c"
#include "t-mul.c"
#include "t-neg.c"
#include "t-sub.c"

test_struct tests[] = {
    TEST_FUNCTION(padic_nmod_add),
    TEST_FUNCTION(padic_nmod_div),
    TEST_FUNCTION(padic_nmod_inv),
    TEST_FUNCTION(padic_nmod_mul),
    TEST_FUNCTION(padic_nmod_neg),
    TEST_FUNCTION(padic_nmod_sub)
};

TEST_MAIN(tests)
