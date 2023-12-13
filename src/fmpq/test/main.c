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

/* For t-get_mpfr.c */
#include <stdio.h>
#include <mpfr.h>

/* Include functions *********************************************************/

#include "t-abs.c"
#include "t-add.c"
#include "t-add_fmpz.c"
#include "t-addmul.c"
#include "t-add_si.c"
#include "t-add_ui.c"
#include "t-canonicalise.c"
#include "t-cfrac_bound.c"
#include "t-cmp.c"
#include "t-dedekind_sum.c"
#include "t-div_2exp.c"
#include "t-div.c"
#include "t-div_fmpz.c"
#include "t-equal_si_ui.c"
#include "t-farey_neighbors.c"
#include "t-gcd_cofactors.c"
#include "t-get_cfrac.c"
#include "t-get_d.c"
#include "t-get_mpfr.c"
#include "t-get_set_str.c"
#include "t-harmonic_ui.c"
#include "t-height.c"
#include "t-init_set_readonly.c"
#include "t-inv.c"
#include "t-mpq_init_set_readonly.c"
#include "t-mul_2exp.c"
#include "t-mul.c"
#include "t-mul_fmpz.c"
#include "t-mul_si.c"
#include "t-mul_ui.c"
#include "t-next_calkin_wilf.c"
#include "t-next_minimal.c"
#include "t-one.c"
#include "t-pow_si.c"
#include "t-randtest.c"
#include "t-reconstruct_fmpz_2.c"
#include "t-reconstruct_fmpz.c"
#include "t-set_cfrac.c"
#include "t-set_fmpz_frac.c"
#include "t-set_si.c"
#include "t-set_ui.c"
#include "t-simplest_between.c"
#include "t-sub.c"
#include "t-sub_fmpz.c"
#include "t-submul.c"
#include "t-sub_si.c"
#include "t-sub_ui.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fmpq_abs),
    TEST_FUNCTION(fmpq_add),
    TEST_FUNCTION(fmpq_add_fmpz),
    TEST_FUNCTION(fmpq_addmul),
    TEST_FUNCTION(fmpq_add_si),
    TEST_FUNCTION(fmpq_add_ui),
    TEST_FUNCTION(fmpq_canonicalise),
    TEST_FUNCTION(fmpq_cfrac_bound),
    TEST_FUNCTION(fmpq_cmp),
    TEST_FUNCTION(fmpq_dedekind_sum),
    TEST_FUNCTION(fmpq_div_2exp),
    TEST_FUNCTION(fmpq_div),
    TEST_FUNCTION(fmpq_div_fmpz),
    TEST_FUNCTION(fmpq_equal_si_ui),
    TEST_FUNCTION(fmpq_farey_neighbors),
    TEST_FUNCTION(fmpq_gcd_cofactors),
    TEST_FUNCTION(fmpq_get_cfrac),
    TEST_FUNCTION(fmpq_get_d),
    TEST_FUNCTION(fmpq_get_mpfr),
    TEST_FUNCTION(fmpq_get_set_str),
    TEST_FUNCTION(fmpq_harmonic_ui),
    TEST_FUNCTION(fmpq_height),
    TEST_FUNCTION(fmpq_init_set_readonly),
    TEST_FUNCTION(fmpq_inv),
    TEST_FUNCTION(fmpq_mpq_init_set_readonly),
    TEST_FUNCTION(fmpq_mul_2exp),
    TEST_FUNCTION(fmpq_mul),
    TEST_FUNCTION(fmpq_mul_fmpz),
    TEST_FUNCTION(fmpq_mul_si),
    TEST_FUNCTION(fmpq_mul_ui),
    TEST_FUNCTION(fmpq_next_calkin_wilf),
    TEST_FUNCTION(fmpq_next_minimal),
    TEST_FUNCTION(fmpq_one),
    TEST_FUNCTION(fmpq_pow_si),
    TEST_FUNCTION(fmpq_randtest),
    TEST_FUNCTION(fmpq_reconstruct_fmpz_2),
    TEST_FUNCTION(fmpq_reconstruct_fmpz),
    TEST_FUNCTION(fmpq_set_cfrac),
    TEST_FUNCTION(fmpq_set_fmpz_frac),
    TEST_FUNCTION(fmpq_set_si),
    TEST_FUNCTION(fmpq_set_ui),
    TEST_FUNCTION(fmpq_simplest_between),
    TEST_FUNCTION(fmpq_sub),
    TEST_FUNCTION(fmpq_sub_fmpz),
    TEST_FUNCTION(fmpq_submul),
    TEST_FUNCTION(fmpq_sub_si),
    TEST_FUNCTION(fmpq_sub_ui)
};

/* main function *************************************************************/

TEST_MAIN(tests)
