/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Try to get fdopen declared for fmpz_mat_[print/read] */
#if defined __STRICT_ANSI__
# undef __STRICT_ANSI__
#endif

#include <string.h>
#include <stdlib.h>

/* Include functions *********************************************************/

#include "t-add_sub.c"
#include "t-can_solve_fflu.c"
#include "t-can_solve_multi_mod_den.c"
#include "t-charpoly_berkowitz.c"
#include "t-charpoly.c"
#include "t-chol_d.c"
#include "t-col_partition.c"
#include "t-concat_horizontal.c"
#include "t-concat_vertical.c"
#include "t-content.c"
#include "t-CRT_ui.c"
#include "t-CRT_ui_unsigned.c"
#include "t-det_bound.c"
#include "t-det.c"
#include "t-det_divisor.c"
#include "t-det_modular_accelerated.c"
#include "t-det_modular.c"
#include "t-entry.c"
#include "t-equal.c"
#include "t-fmpz_vec_mul.c"
#include "t-get_d_mat.c"
#include "t-get_d_mat_transpose.c"
#include "t-get_nmod_mat.c"
#include "t-gram.c"
#include "t-hadamard.c"
#include "t-hnf.c"
#include "t-hnf_classical.c"
#include "t-hnf_minors.c"
#include "t-hnf_minors_transform.c"
#include "t-hnf_modular.c"
#include "t-hnf_modular_eldiv.c"
#include "t-hnf_pernet_stein.c"
#include "t-hnf_transform.c"
#include "t-hnf_xgcd.c"
#include "t-howell_form_mod.c"
#include "t-init_clear.c"
#include "t-inv.c"
#include "t-invert_rows_cols.c"
#include "t-is_empty.c"
#include "t-is_one.c"
#include "t-is_spd.c"
#include "t-is_square.c"
#include "t-is_zero.c"
#include "t-kronecker_product.c"
#include "t-lll_original.c"
#include "t-lll_storjohann.c"
#include "t-max_bits.c"
#include "t-minpoly.c"
#include "t-mul_blas.c"
#include "t-mul.c"
#include "t-mul_classical.c"
#include "t-mul_double_word.c"
#include "t-mul_fft.c"
#include "t-mul_fmpz_vec.c"
#include "t-mul_multi_mod.c"
#include "t-mul_small.c"
#include "t-mul_strassen.c"
#include "t-multi_CRT_ui.c"
#include "t-multi_CRT_ui_unsigned.c"
#include "t-nullspace.c"
#include "t-one.c"
#include "t-pow.c"
#include "t-print_read.c"
#include "t-rank.c"
#include "t-rref.c"
#include "t-rref_fflu.c"
#include "t-rref_mod.c"
#include "t-rref_mul.c"
#include "t-scalar_addmul_fmpz.c"
#include "t-scalar_addmul_nmod_mat_fmpz.c"
#include "t-scalar_addmul_nmod_mat_ui.c"
#include "t-scalar_addmul_si.c"
#include "t-scalar_addmul_ui.c"
#include "t-scalar_mod_fmpz.c"
#include "t-scalar_mul_2exp.c"
#include "t-scalar_mul_fmpz.c"
#include "t-scalar_mul_si.c"
#include "t-scalar_mul_ui.c"
#include "t-scalar_smod.c"
#include "t-snf_diagonal.c"
#include "t-snf_iliopoulos.c"
#include "t-snf_kannan_bachem.c"
#include "t-solve_bound.c"
#include "t-solve.c"
#include "t-solve_cramer.c"
#include "t-solve_dixon.c"
#include "t-solve_dixon_den.c"
#include "t-solve_fflu.c"
#include "t-solve_multi_mod_den.c"
#include "t-sqr.c"
#include "t-trace.c"
#include "t-transpose.c"
#include "t-window_init_clear.c"
#include "t-zero.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(fmpz_mat_add_sub),
    TEST_FUNCTION(fmpz_mat_can_solve_fflu),
    TEST_FUNCTION(fmpz_mat_can_solve_multi_mod_den),
    TEST_FUNCTION(fmpz_mat_charpoly_berkowitz),
    TEST_FUNCTION(fmpz_mat_charpoly),
    TEST_FUNCTION(fmpz_mat_chol_d),
    TEST_FUNCTION(fmpz_mat_col_partition),
    TEST_FUNCTION(fmpz_mat_concat_horizontal),
    TEST_FUNCTION(fmpz_mat_concat_vertical),
    TEST_FUNCTION(fmpz_mat_content),
    TEST_FUNCTION(fmpz_mat_CRT_ui),
    TEST_FUNCTION(fmpz_mat_CRT_ui_unsigned),
    TEST_FUNCTION(fmpz_mat_det_bound),
    TEST_FUNCTION(fmpz_mat_det),
    TEST_FUNCTION(fmpz_mat_det_divisor),
    TEST_FUNCTION(fmpz_mat_det_modular_accelerated),
    TEST_FUNCTION(fmpz_mat_det_modular),
    TEST_FUNCTION(fmpz_mat_entry),
    TEST_FUNCTION(fmpz_mat_equal),
    TEST_FUNCTION(fmpz_mat_fmpz_vec_mul),
    TEST_FUNCTION(fmpz_mat_get_d_mat),
    TEST_FUNCTION(fmpz_mat_get_d_mat_transpose),
    TEST_FUNCTION(fmpz_mat_get_nmod_mat),
    TEST_FUNCTION(fmpz_mat_gram),
    TEST_FUNCTION(fmpz_mat_hadamard),
    TEST_FUNCTION(fmpz_mat_hnf),
    TEST_FUNCTION(fmpz_mat_hnf_classical),
    TEST_FUNCTION(fmpz_mat_hnf_minors),
    TEST_FUNCTION(fmpz_mat_hnf_minors_transform),
    TEST_FUNCTION(fmpz_mat_hnf_modular),
    TEST_FUNCTION(fmpz_mat_hnf_modular_eldiv),
    TEST_FUNCTION(fmpz_mat_hnf_pernet_stein),
    TEST_FUNCTION(fmpz_mat_hnf_transform),
    TEST_FUNCTION(fmpz_mat_hnf_xgcd),
    TEST_FUNCTION(fmpz_mat_howell_form_mod),
    TEST_FUNCTION(fmpz_mat_init_clear),
    TEST_FUNCTION(fmpz_mat_inv),
    TEST_FUNCTION(fmpz_mat_invert_rows_cols),
    TEST_FUNCTION(fmpz_mat_is_empty),
    TEST_FUNCTION(fmpz_mat_is_one),
    TEST_FUNCTION(fmpz_mat_is_spd),
    TEST_FUNCTION(fmpz_mat_is_square),
    TEST_FUNCTION(fmpz_mat_is_zero),
    TEST_FUNCTION(fmpz_mat_kronecker_product),
    TEST_FUNCTION(fmpz_mat_lll_original),
    TEST_FUNCTION(fmpz_mat_lll_storjohann),
    TEST_FUNCTION(fmpz_mat_max_bits),
    TEST_FUNCTION(fmpz_mat_minpoly),
    TEST_FUNCTION(fmpz_mat_mul_blas),
    TEST_FUNCTION(fmpz_mat_mul),
    TEST_FUNCTION(fmpz_mat_mul_classical),
    TEST_FUNCTION(fmpz_mat_mul_double_word),
    TEST_FUNCTION(fmpz_mat_mul_fft),
    TEST_FUNCTION(fmpz_mat_mul_fmpz_vec),
    TEST_FUNCTION(fmpz_mat_mul_multi_mod),
    TEST_FUNCTION(fmpz_mat_mul_small),
    TEST_FUNCTION(fmpz_mat_mul_strassen),
    TEST_FUNCTION(fmpz_mat_multi_CRT_ui),
    TEST_FUNCTION(fmpz_mat_multi_CRT_ui_unsigned),
    TEST_FUNCTION(fmpz_mat_nullspace),
    TEST_FUNCTION(fmpz_mat_one),
    TEST_FUNCTION(fmpz_mat_pow),
    TEST_FUNCTION(fmpz_mat_print_read),
    TEST_FUNCTION(fmpz_mat_rank),
    TEST_FUNCTION(fmpz_mat_rref),
    TEST_FUNCTION(fmpz_mat_rref_fflu),
    TEST_FUNCTION(fmpz_mat_rref_mod),
    TEST_FUNCTION(fmpz_mat_rref_mul),
    TEST_FUNCTION(fmpz_mat_scalar_addmul_fmpz),
    TEST_FUNCTION(fmpz_mat_scalar_addmul_nmod_mat_fmpz),
    TEST_FUNCTION(fmpz_mat_scalar_addmul_nmod_mat_ui),
    TEST_FUNCTION(fmpz_mat_scalar_addmul_si),
    TEST_FUNCTION(fmpz_mat_scalar_addmul_ui),
    TEST_FUNCTION(fmpz_mat_scalar_mod_fmpz),
    TEST_FUNCTION(fmpz_mat_scalar_mul_2exp),
    TEST_FUNCTION(fmpz_mat_scalar_mul_fmpz),
    TEST_FUNCTION(fmpz_mat_scalar_mul_si),
    TEST_FUNCTION(fmpz_mat_scalar_mul_ui),
    TEST_FUNCTION(fmpz_mat_scalar_smod),
    TEST_FUNCTION(fmpz_mat_snf_diagonal),
    TEST_FUNCTION(fmpz_mat_snf_iliopoulos),
    TEST_FUNCTION(fmpz_mat_snf_kannan_bachem),
    TEST_FUNCTION(fmpz_mat_solve_bound),
    TEST_FUNCTION(fmpz_mat_solve),
    TEST_FUNCTION(fmpz_mat_solve_cramer),
    TEST_FUNCTION(fmpz_mat_solve_dixon),
    TEST_FUNCTION(fmpz_mat_solve_dixon_den),
    TEST_FUNCTION(fmpz_mat_solve_fflu),
    TEST_FUNCTION(fmpz_mat_solve_multi_mod_den),
    TEST_FUNCTION(fmpz_mat_sqr),
    TEST_FUNCTION(fmpz_mat_trace),
    TEST_FUNCTION(fmpz_mat_transpose),
    TEST_FUNCTION(fmpz_mat_window_init_clear),
    TEST_FUNCTION(fmpz_mat_zero)
};

/* main function *************************************************************/

TEST_MAIN(tests)
