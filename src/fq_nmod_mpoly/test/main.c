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

#include "t-add_sub.c"
#include "t-add_sub_fq_nmod.c"
#include "t-cmp.c"
#include "t-compose_fq_nmod_mpoly.c"
#include "t-compose_fq_nmod_poly.c"
#include "t-degree.c"
#include "t-derivative.c"
#include "t-div_monagan_pearce.c"
#include "t-divrem_ideal_monagan_pearce.c"
#include "t-divrem_monagan_pearce.c"
#include "t-evaluate.c"
#include "t-gcd_brown.c"
#include "t-gcd.c"
#include "t-gcd_cofactors.c"
#include "t-gcd_hensel.c"
#include "t-gcd_zippel2.c"
#include "t-gcd_zippel.c"
#include "t-gen.c"
#include "t-get_coeff_vars_ui.c"
#include "t-get_set_coeff.c"
#include "t-get_set_coeff_fq_nmod_monomial.c"
#include "t-get_set_is_fq_nmod.c"
#include "t-get_set_str_pretty.c"
#include "t-get_set_term_coeff_fq_nmod.c"
#include "t-get_set_term_exp_fmpz.c"
#include "t-get_set_term_exp_si.c"
#include "t-get_set_term_exp_ui.c"
#include "t-get_term.c"
#include "t-get_term_monomial.c"
#include "t-mpolyuu_divides.c"
#include "t-mul_johnson.c"
#include "t-push_term_fq_nmod_fmpz.c"
#include "t-push_term_fq_nmod_ui.c"
#include "t-quadratic_root.c"
#include "t-repack_bits.c"
#include "t-resize.c"
#include "t-resultant_discriminant.c"
#include "t-reverse.c"
#include "t-scalar_addmul_fq_nmod.c"
#include "t-scalar_mul_fq_nmod.c"
#include "t-sort_terms.c"
#include "t-sqrt.c"
#include "t-total_degree.c"
#include "t-univar.c"
#include "t-univar_resultant.c"
#include "t-used_vars.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fq_nmod_mpoly_add_sub),
    TEST_FUNCTION(fq_nmod_mpoly_add_sub_fq_nmod),
    TEST_FUNCTION(fq_nmod_mpoly_cmp),
    TEST_FUNCTION(fq_nmod_mpoly_compose_fq_nmod_mpoly),
    TEST_FUNCTION(fq_nmod_mpoly_compose_fq_nmod_poly),
    TEST_FUNCTION(fq_nmod_mpoly_degree),
    TEST_FUNCTION(fq_nmod_mpoly_derivative),
    TEST_FUNCTION(fq_nmod_mpoly_div_monagan_pearce),
    TEST_FUNCTION(fq_nmod_mpoly_divrem_ideal_monagan_pearce),
    TEST_FUNCTION(fq_nmod_mpoly_divrem_monagan_pearce),
    TEST_FUNCTION(fq_nmod_mpoly_evaluate),
    TEST_FUNCTION(fq_nmod_mpoly_gcd_brown),
    TEST_FUNCTION(fq_nmod_mpoly_gcd),
    TEST_FUNCTION(fq_nmod_mpoly_gcd_cofactors),
    TEST_FUNCTION(fq_nmod_mpoly_gcd_hensel),
    TEST_FUNCTION(fq_nmod_mpoly_gcd_zippel2),
    TEST_FUNCTION(fq_nmod_mpoly_gcd_zippel),
    TEST_FUNCTION(fq_nmod_mpoly_gen),
    TEST_FUNCTION(fq_nmod_mpoly_get_coeff_vars_ui),
    TEST_FUNCTION(fq_nmod_mpoly_get_set_coeff),
    TEST_FUNCTION(fq_nmod_mpoly_get_set_coeff_fq_nmod_monomial),
    TEST_FUNCTION(fq_nmod_mpoly_get_set_is_fq_nmod),
    TEST_FUNCTION(fq_nmod_mpoly_get_set_str_pretty),
    TEST_FUNCTION(fq_nmod_mpoly_get_set_term_coeff_fq_nmod),
    TEST_FUNCTION(fq_nmod_mpoly_get_set_term_exp_fmpz),
    TEST_FUNCTION(fq_nmod_mpoly_get_set_term_exp_si),
    TEST_FUNCTION(fq_nmod_mpoly_get_set_term_exp_ui),
    TEST_FUNCTION(fq_nmod_mpoly_get_term),
    TEST_FUNCTION(fq_nmod_mpoly_get_term_monomial),
    TEST_FUNCTION(fq_nmod_mpoly_mpolyuu_divides),
    TEST_FUNCTION(fq_nmod_mpoly_mul_johnson),
    TEST_FUNCTION(fq_nmod_mpoly_push_term_fq_nmod_fmpz),
    TEST_FUNCTION(fq_nmod_mpoly_push_term_fq_nmod_ui),
    TEST_FUNCTION(fq_nmod_mpoly_quadratic_root),
    TEST_FUNCTION(fq_nmod_mpoly_repack_bits),
    TEST_FUNCTION(fq_nmod_mpoly_resize),
    TEST_FUNCTION(fq_nmod_mpoly_resultant_discriminant),
    TEST_FUNCTION(fq_nmod_mpoly_reverse),
    TEST_FUNCTION(fq_nmod_mpoly_scalar_addmul_fq_nmod),
    TEST_FUNCTION(fq_nmod_mpoly_scalar_mul_fq_nmod),
    TEST_FUNCTION(fq_nmod_mpoly_sort_terms),
    TEST_FUNCTION(fq_nmod_mpoly_sqrt),
    TEST_FUNCTION(fq_nmod_mpoly_total_degree),
    TEST_FUNCTION(fq_nmod_mpoly_univar),
    TEST_FUNCTION(fq_nmod_mpoly_univar_resultant),
    TEST_FUNCTION(fq_nmod_mpoly_used_vars)
};

/* main function *************************************************************/

TEST_MAIN(tests)
