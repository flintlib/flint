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

#include "mpoly.h"
#include "nf.h"
#include "gr.h"

/* Include functions *********************************************************/

#include "t-acb.c"
#include "t-arb.c"
#include "t-ca.c"
#include "t-dirichlet.c"
#include "t-fmpq.c"
#include "t-fmpq_poly.c"
#include "t-fmpz.c"
#include "t-fmpzi.c"
#include "t-fmpz_mod.c"
#include "t-fmpz_mpoly.c"
#include "t-fmpz_mpoly_q.c"
#include "t-fmpz_poly.c"
#include "t-fq.c"
#include "t-fq_nmod.c"
#include "t-fq_zech.c"
#include "t-matrix_acb.c"
#include "t-matrix_arb.c"
#include "t-matrix_fmpq.c"
#include "t-matrix_fmpz.c"
#include "t-matrix_nmod8.c"
#include "t-mpoly_nmod8.c"
#include "t-nf.c"
#include "t-nmod32.c"
#include "t-nmod8.c"
#include "t-nmod.c"
#include "t-perm.c"
#include "t-polynomial_acb.c"
#include "t-polynomial_arb.c"
#include "t-polynomial_fmpq.c"
#include "t-polynomial_fmpz.c"
#include "t-polynomial_nmod8.c"
#include "t-psl2z.c"
#include "t-qqbar.c"
#include "t-series_acb.c"
#include "t-series_arb.c"
#include "t-series_fmpq.c"
#include "t-series_fmpz.c"
#include "t-series_nmod8.c"
#include "t-vector_acb.c"
#include "t-vector_arb.c"
#include "t-vector_fmpz.c"
#include "t-vector_fmpq.c"
#include "t-vector_nmod.c"
#include "t-vector_nmod8.c"
#include "t-vector_nmod32.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(gr_acb),
    TEST_FUNCTION(gr_arb),
    TEST_FUNCTION(gr_ca),
    TEST_FUNCTION(gr_dirichlet),
    TEST_FUNCTION(gr_fmpq),
    TEST_FUNCTION(gr_fmpq_poly),
    TEST_FUNCTION(gr_fmpz),
    TEST_FUNCTION(gr_fmpzi),
    TEST_FUNCTION(gr_fmpz_mod),
    TEST_FUNCTION(gr_fmpz_mpoly),
    TEST_FUNCTION(gr_fmpz_mpoly_q),
    TEST_FUNCTION(gr_fmpz_poly),
    TEST_FUNCTION(gr_fq),
    TEST_FUNCTION(gr_fq_nmod),
    TEST_FUNCTION(gr_fq_zech),
    TEST_FUNCTION(gr_matrix_acb),
    TEST_FUNCTION(gr_matrix_arb),
    TEST_FUNCTION(gr_matrix_fmpq),
    TEST_FUNCTION(gr_matrix_fmpz),
    TEST_FUNCTION(gr_matrix_nmod8),
    TEST_FUNCTION(gr_mpoly_nmod8),
    TEST_FUNCTION(gr_nf),
    TEST_FUNCTION(gr_nmod32),
    TEST_FUNCTION(gr_nmod8),
    TEST_FUNCTION(gr_nmod),
    TEST_FUNCTION(gr_perm),
    TEST_FUNCTION(gr_polynomial_acb),
    TEST_FUNCTION(gr_polynomial_arb),
    TEST_FUNCTION(gr_polynomial_fmpq),
    TEST_FUNCTION(gr_polynomial_fmpz),
    TEST_FUNCTION(gr_polynomial_nmod8),
    TEST_FUNCTION(gr_psl2z),
    TEST_FUNCTION(gr_qqbar),
    TEST_FUNCTION(gr_series_acb),
    TEST_FUNCTION(gr_series_arb),
    TEST_FUNCTION(gr_series_fmpq),
    TEST_FUNCTION(gr_series_fmpz),
    TEST_FUNCTION(gr_series_nmod8),
    TEST_FUNCTION(gr_vector_acb),
    TEST_FUNCTION(gr_vector_arb),
    TEST_FUNCTION(gr_vector_fmpz),
    TEST_FUNCTION(gr_vector_fmpq),
    TEST_FUNCTION(gr_vector_nmod),
    TEST_FUNCTION(gr_vector_nmod8),
    TEST_FUNCTION(gr_vector_nmod32),
};

/* main function *************************************************************/

TEST_MAIN(tests)
