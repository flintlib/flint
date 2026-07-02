/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-add_sub.c"
#include "t-gen.c"
#include "t-get_set_coeff.c"
#include "t-integral.c"
#include "t-divides.c"
#include "t-divrem.c"
#include "t-divrem_ideal.c"
#include "t-quasidivrem.c"
#include "t-quasidivrem_ideal.c"
#include "t-vec.c"
#include "t-mul_heap.c"
#include "t-mul_heap_threaded.c"
#include "t-mul_monomial.c"
#include "t-sqr.c"
#include "t-sqr_commutative_heap.c"
#include "t-sqr_commutative_heap_threaded.c"
#include "t-ring.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_mpoly_add_sub),
    TEST_FUNCTION(gr_mpoly_gen),
    TEST_FUNCTION(gr_mpoly_get_set_coeff),
    TEST_FUNCTION(gr_mpoly_integral),
    TEST_FUNCTION(gr_mpoly_divides),
    TEST_FUNCTION(gr_mpoly_divrem),
    TEST_FUNCTION(gr_mpoly_divrem_ideal),
    TEST_FUNCTION(gr_mpoly_quasidivrem),
    TEST_FUNCTION(gr_mpoly_quasidivrem_ideal),
    TEST_FUNCTION(gr_mpoly_vec),
    TEST_FUNCTION(gr_mpoly_mul_heap),
    TEST_FUNCTION(gr_mpoly_mul_heap_threaded),
    TEST_FUNCTION(gr_mpoly_mul_monomial),
    TEST_FUNCTION(gr_mpoly_ring),
    TEST_FUNCTION(gr_mpoly_sqr),
    TEST_FUNCTION(gr_mpoly_sqr_commutative_heap),
    TEST_FUNCTION(gr_mpoly_sqr_commutative_heap_threaded),
};

/* main function *************************************************************/

TEST_MAIN(tests)
