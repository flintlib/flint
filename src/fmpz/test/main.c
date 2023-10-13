/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Try to get fdopen declared for fmpz_[print/read] */
#if defined __STRICT_ANSI__
# undef __STRICT_ANSI__
#endif

#include <string.h>
#include <stdio.h>
#include <mpfr.h>

/* Include functions *********************************************************/

#include "t-abs.c"
#include "t-abs_fits_ui.c"
#include "t-abs_lbound_ui_2exp.c"
#include "t-abs_ubound_ui_2exp.c"
#include "t-add.c"
#include "t-addmul.c"
#include "t-addmul_si.c"
#include "t-addmul_ui.c"
#include "t-add_ui.c"
#include "t-and.c"
#include "t-bin_uiui.c"
#include "t-bit_pack.c"
#include "t-bits.c"
#include "t-cdiv_q_2exp.c"
#include "t-cdiv_q.c"
#include "t-cdiv_qr.c"
#include "t-cdiv_q_si.c"
#include "t-cdiv_q_ui.c"
#include "t-cdiv_r_2exp.c"
#include "t-cdiv_ui.c"
#include "t-clog.c"
#include "t-clog_ui.c"
#include "t-cmp2abs.c"
#include "t-cmpabs.c"
#include "t-cmp.c"
#include "t-cmp_si.c"
#include "t-cmp_ui.c"
#include "t-comb_init_clear.c"
#include "t-combit.c"
#include "t-complement.c"
#include "t-crt.c"
#include "t-crt_ui.c"
#include "t-divexact2_uiui.c"
#include "t-divexact.c"
#include "t-divexact_si.c"
#include "t-divexact_ui.c"
#include "t-divides.c"
#include "t-divides_mod_list.c"
#include "t-divisible.c"
#include "t-divisible_si.c"
#include "t-divisor_in_residue_class_lenstra.c"
#include "t-divisor_sigma.c"
#include "t-div_newton.c"
#include "t-dlog.c"
#include "t-equal.c"
#include "t-equal_si.c"
#include "t-equal_ui.c"
#include "t-euler_phi.c"
#include "t-fac_ui.c"
#include "t-fdiv_q_2exp.c"
#include "t-fdiv_q.c"
#include "t-fdiv_qr.c"
#include "t-fdiv_qr_preinvn.c"
#include "t-fdiv_q_si.c"
#include "t-fdiv_q_ui.c"
#include "t-fdiv_r_2exp.c"
#include "t-fdiv_r.c"
#include "t-fdiv_ui.c"
#include "t-fib_ui.c"
#include "t-fits_si.c"
#include "t-flog.c"
#include "t-flog_ui.c"
#include "t-fmma.c"
#include "t-fmms.c"
#include "t-fmpz.c"
#include "t-fmpz_cleanup.c"
#include "t-fmpz_stress.c"
#include "t-gcd3.c"
#include "t-gcd.c"
#include "t-gcdinv.c"
#include "t-gcd_ui.c"
#include "t-get_d_2exp.c"
#include "t-get_d.c"
#include "t-get_mpf.c"
#include "t-get_mpfr.c"
#include "t-get_mpn.c"
#include "t-get_mpz.c"
#include "t-get_nmod.c"
#include "t-get_set_ui_array.c"
#include "t-get_si.c"
#include "t-get_str.c"
#include "t-get_ui.c"
#include "t-init2.c"
#include "t-init_set.c"
#include "t-init_set_readonly.c"
#include "t-init_set_ui.c"
#include "t-invmod.c"
#include "t-is_even.c"
#include "t-is_perfect_power.c"
#include "t-is_prime.c"
#include "t-is_prime_morrison.c"
#include "t-is_prime_pocklington.c"
#include "t-is_prime_pseudosquare.c"
#include "t-is_probabprime_BPSW.c"
#include "t-is_probabprime_lucas.c"
#include "t-is_square.c"
#include "t-is_strong_probabprime.c"
#include "t-jacobi.c"
#include "t-kronecker.c"
#include "t-lcm.c"
#include "t-mod.c"
#include "t-mod_ui.c"
#include "t-moebius_mu.c"
#include "t-mpz_init_set_readonly.c"
#include "t-mul_2exp.c"
#include "t-mul2_uiui.c"
#include "t-mul.c"
#include "t-mul_si.c"
#include "t-mul_si_tdiv_q_2exp.c"
#include "t-mul_tdiv_q_2exp.c"
#include "t-multi_CRT_multi_mod.c"
#include "t-multi_CRT_ui.c"
#include "t-mul_ui.c"
#include "t-ndiv_qr.c"
#include "t-neg.c"
#include "t-neg_ui.c"
#include "t-neg_uiui.c"
#include "t-nextprime.c"
#include "t-or.c"
#include "t-out_inp_raw.c"
#include "t-popcnt.c"
#include "t-powm.c"
#include "t-powm_ui.c"
#include "t-pow_ui.c"
#include "t-primorial.c"
#include "t-print_read.c"
#include "t-randprime.c"
#include "t-remove.c"
#include "t-rfac_ui.c"
#include "t-rfac_uiui.c"
#include "t-root.c"
#include "t-setbit.c"
#include "t-set.c"
#include "t-set_d_2exp.c"
#include "t-set_signed_ui_array.c"
#include "t-set_signed_uiui.c"
#include "t-set_signed_uiuiui.c"
#include "t-set_str.c"
#include "t-set_ui_smod.c"
#include "t-set_uiui.c"
#include "t-sgn.c"
#include "t-size.c"
#include "t-sizeinbase.c"
#include "t-smod.c"
#include "t-sqrt.c"
#include "t-sqrtmod.c"
#include "t-sqrtrem.c"
#include "t-sub.c"
#include "t-submul.c"
#include "t-submul_si.c"
#include "t-submul_ui.c"
#include "t-sub_ui.c"
#include "t-swap.c"
#include "t-tdiv_q_2exp.c"
#include "t-tdiv_q.c"
#include "t-tdiv_qr.c"
#include "t-tdiv_q_si.c"
#include "t-tdiv_q_ui.c"
#include "t-tdiv_r_2exp.c"
#include "t-tdiv_ui.c"
#include "t-tstbit.c"
#include "t-val2.c"
#include "t-xgcd.c"
#include "t-xgcd_canonical_bezout.c"
#include "t-xgcd_partial.c"
#include "t-xor.c"

/* Array of test functions ***************************************************/

int (*test_functions[])(void) =
{
    TEST_FUNCTION(fmpz_abs),
    TEST_FUNCTION(fmpz_abs_fits_ui),
    TEST_FUNCTION(fmpz_abs_lbound_ui_2exp),
    TEST_FUNCTION(fmpz_abs_ubound_ui_2exp),
    TEST_FUNCTION(fmpz_add),
    TEST_FUNCTION(fmpz_addmul),
    TEST_FUNCTION(fmpz_addmul_si),
    TEST_FUNCTION(fmpz_addmul_ui),
    TEST_FUNCTION(fmpz_add_ui),
    TEST_FUNCTION(fmpz_and),
    TEST_FUNCTION(fmpz_bin_uiui),
    TEST_FUNCTION(fmpz_bit_pack),
    TEST_FUNCTION(fmpz_bits),
    TEST_FUNCTION(fmpz_cdiv_q_2exp),
    TEST_FUNCTION(fmpz_cdiv_q),
    TEST_FUNCTION(fmpz_cdiv_qr),
    TEST_FUNCTION(fmpz_cdiv_q_si),
    TEST_FUNCTION(fmpz_cdiv_q_ui),
    TEST_FUNCTION(fmpz_cdiv_r_2exp),
    TEST_FUNCTION(fmpz_cdiv_ui),
    TEST_FUNCTION(fmpz_clog),
    TEST_FUNCTION(fmpz_clog_ui),
    TEST_FUNCTION(fmpz_cmp2abs),
    TEST_FUNCTION(fmpz_cmpabs),
    TEST_FUNCTION(fmpz_cmp),
    TEST_FUNCTION(fmpz_cmp_si),
    TEST_FUNCTION(fmpz_cmp_ui),
    TEST_FUNCTION(fmpz_comb_init_clear),
    TEST_FUNCTION(fmpz_combit),
    TEST_FUNCTION(fmpz_complement),
    TEST_FUNCTION(fmpz_CRT),
    TEST_FUNCTION(fmpz_CRT_ui),
    TEST_FUNCTION(fmpz_divexact2_uiui),
    TEST_FUNCTION(fmpz_divexact),
    TEST_FUNCTION(fmpz_divexact_si),
    TEST_FUNCTION(fmpz_divexact_ui),
    TEST_FUNCTION(fmpz_divides),
    TEST_FUNCTION(fmpz_divides_mod_list),
    TEST_FUNCTION(fmpz_divisible),
    TEST_FUNCTION(fmpz_divisible_si),
    TEST_FUNCTION(fmpz_divisor_in_residue_class_lenstra),
    TEST_FUNCTION(fmpz_divisor_sigma),
    TEST_FUNCTION(fmpz_div_newton),
    TEST_FUNCTION(fmpz_dlog),
    TEST_FUNCTION(fmpz_equal),
    TEST_FUNCTION(fmpz_equal_si),
    TEST_FUNCTION(fmpz_equal_ui),
    TEST_FUNCTION(fmpz_euler_phi),
    TEST_FUNCTION(fmpz_fac_ui),
    TEST_FUNCTION(fmpz_fdiv_q_2exp),
    TEST_FUNCTION(fmpz_fdiv_q),
    TEST_FUNCTION(fmpz_fdiv_qr),
    TEST_FUNCTION(fmpz_fdiv_qr_preinvn),
    TEST_FUNCTION(fmpz_fdiv_q_si),
    TEST_FUNCTION(fmpz_fdiv_q_ui),
    TEST_FUNCTION(fmpz_fdiv_r_2exp),
    TEST_FUNCTION(fmpz_fdiv_r),
    TEST_FUNCTION(fmpz_fdiv_ui),
    TEST_FUNCTION(fmpz_fib_ui),
    TEST_FUNCTION(fmpz_fits_si),
    TEST_FUNCTION(fmpz_flog),
    TEST_FUNCTION(fmpz_flog_ui),
    TEST_FUNCTION(fmpz_fmma),
    TEST_FUNCTION(fmpz_fmms),
    TEST_FUNCTION(fmpz_fmpz),
    TEST_FUNCTION(fmpz_cleanup),
    TEST_FUNCTION(fmpz_stress),
    TEST_FUNCTION(fmpz_gcd3),
    TEST_FUNCTION(fmpz_gcd),
    TEST_FUNCTION(fmpz_gcdinv),
    TEST_FUNCTION(fmpz_gcd_ui),
    TEST_FUNCTION(fmpz_get_d_2exp),
    TEST_FUNCTION(fmpz_get_d),
    TEST_FUNCTION(fmpz_get_mpf),
    TEST_FUNCTION(fmpz_get_mpfr),
    TEST_FUNCTION(fmpz_get_mpn),
    TEST_FUNCTION(fmpz_get_mpz),
    TEST_FUNCTION(fmpz_get_nmod),
    TEST_FUNCTION(fmpz_get_set_ui_array),
    TEST_FUNCTION(fmpz_get_si),
    TEST_FUNCTION(fmpz_get_str),
    TEST_FUNCTION(fmpz_get_ui),
    TEST_FUNCTION(fmpz_init2),
    TEST_FUNCTION(fmpz_init_set),
    TEST_FUNCTION(fmpz_init_set_readonly),
    TEST_FUNCTION(fmpz_init_set_ui),
    TEST_FUNCTION(fmpz_invmod),
    TEST_FUNCTION(fmpz_is_even),
    TEST_FUNCTION(fmpz_is_perfect_power),
    TEST_FUNCTION(fmpz_is_prime),
    TEST_FUNCTION(fmpz_is_prime_morrison),
    TEST_FUNCTION(fmpz_is_prime_pocklington),
    TEST_FUNCTION(fmpz_is_prime_pseudosquare),
    TEST_FUNCTION(fmpz_is_probabprime_BPSW),
    TEST_FUNCTION(fmpz_is_probabprime_lucas),
    TEST_FUNCTION(fmpz_is_square),
    TEST_FUNCTION(fmpz_is_strong_probabprime),
    TEST_FUNCTION(fmpz_jacobi),
    TEST_FUNCTION(fmpz_kronecker),
    TEST_FUNCTION(fmpz_lcm),
    TEST_FUNCTION(fmpz_mod),
    TEST_FUNCTION(fmpz_mod_ui),
    TEST_FUNCTION(fmpz_moebius_mu),
    TEST_FUNCTION(fmpz_mpz_init_set_readonly),
    TEST_FUNCTION(fmpz_mul_2exp),
    TEST_FUNCTION(fmpz_mul2_uiui),
    TEST_FUNCTION(fmpz_mul),
    TEST_FUNCTION(fmpz_mul_si),
    TEST_FUNCTION(fmpz_mul_si_tdiv_q_2exp),
    TEST_FUNCTION(fmpz_mul_tdiv_q_2exp),
    TEST_FUNCTION(fmpz_multi_CRT_multi_mod),
    TEST_FUNCTION(fmpz_multi_CRT_ui),
    TEST_FUNCTION(fmpz_mul_ui),
    TEST_FUNCTION(fmpz_ndiv_qr),
    TEST_FUNCTION(fmpz_neg),
    TEST_FUNCTION(fmpz_neg_ui),
    TEST_FUNCTION(fmpz_neg_uiui),
    TEST_FUNCTION(fmpz_nextprime),
    TEST_FUNCTION(fmpz_or),
    TEST_FUNCTION(fmpz_out_inp_raw),
    TEST_FUNCTION(fmpz_popcnt),
    TEST_FUNCTION(fmpz_powm),
    TEST_FUNCTION(fmpz_powm_ui),
    TEST_FUNCTION(fmpz_pow_ui),
    TEST_FUNCTION(fmpz_primorial),
    TEST_FUNCTION(fmpz_print_read),
    TEST_FUNCTION(fmpz_randprime),
    TEST_FUNCTION(fmpz_remove),
    TEST_FUNCTION(fmpz_rfac_ui),
    TEST_FUNCTION(fmpz_rfac_uiui),
    TEST_FUNCTION(fmpz_root),
    TEST_FUNCTION(fmpz_setbit),
    TEST_FUNCTION(fmpz_set),
    TEST_FUNCTION(fmpz_set_d_2exp),
    TEST_FUNCTION(fmpz_set_signed_ui_array),
    TEST_FUNCTION(fmpz_set_signed_uiui),
    TEST_FUNCTION(fmpz_set_signed_uiuiui),
    TEST_FUNCTION(fmpz_set_str),
    TEST_FUNCTION(fmpz_set_ui_smod),
    TEST_FUNCTION(fmpz_set_uiui),
    TEST_FUNCTION(fmpz_sgn),
    TEST_FUNCTION(fmpz_size),
    TEST_FUNCTION(fmpz_sizeinbase),
    TEST_FUNCTION(fmpz_smod),
    TEST_FUNCTION(fmpz_sqrt),
    TEST_FUNCTION(fmpz_sqrtmod),
    TEST_FUNCTION(fmpz_sqrtrem),
    TEST_FUNCTION(fmpz_sub),
    TEST_FUNCTION(fmpz_submul),
    TEST_FUNCTION(fmpz_submul_si),
    TEST_FUNCTION(fmpz_submul_ui),
    TEST_FUNCTION(fmpz_sub_ui),
    TEST_FUNCTION(fmpz_swap),
    TEST_FUNCTION(fmpz_tdiv_q_2exp),
    TEST_FUNCTION(fmpz_tdiv_q),
    TEST_FUNCTION(fmpz_tdiv_qr),
    TEST_FUNCTION(fmpz_tdiv_q_si),
    TEST_FUNCTION(fmpz_tdiv_q_ui),
    TEST_FUNCTION(fmpz_tdiv_r_2exp),
    TEST_FUNCTION(fmpz_tdiv_ui),
    TEST_FUNCTION(fmpz_tstbit),
    TEST_FUNCTION(fmpz_val2),
    TEST_FUNCTION(fmpz_xgcd),
    TEST_FUNCTION(fmpz_xgcd_canonical_bezout),
    TEST_FUNCTION(fmpz_xgcd_partial),
    TEST_FUNCTION(fmpz_xor)
};

char fmpz_abs_name[] = "fmpz_abs";
char fmpz_abs_fits_ui_name[] = "fmpz_abs_fits_ui";
char fmpz_abs_lbound_ui_2exp_name[] = "fmpz_abs_lbound_ui_2exp";
char fmpz_abs_ubound_ui_2exp_name[] = "fmpz_abs_ubound_ui_2exp";
char fmpz_add_name[] = "fmpz_add";
char fmpz_addmul_name[] = "fmpz_addmul";
char fmpz_addmul_si_name[] = "fmpz_addmul_si";
char fmpz_addmul_ui_name[] = "fmpz_addmul_ui";
char fmpz_add_ui_name[] = "fmpz_add_ui";
char fmpz_and_name[] = "fmpz_and";
char fmpz_bin_uiui_name[] = "fmpz_bin_uiui";
char fmpz_bit_pack_name[] = "fmpz_bit_pack";
char fmpz_bits_name[] = "fmpz_bits";
char fmpz_cdiv_q_2exp_name[] = "fmpz_cdiv_q_2exp";
char fmpz_cdiv_q_name[] = "fmpz_cdiv_q";
char fmpz_cdiv_qr_name[] = "fmpz_cdiv_qr";
char fmpz_cdiv_q_si_name[] = "fmpz_cdiv_q_si";
char fmpz_cdiv_q_ui_name[] = "fmpz_cdiv_q_ui";
char fmpz_cdiv_r_2exp_name[] = "fmpz_cdiv_r_2exp";
char fmpz_cdiv_ui_name[] = "fmpz_cdiv_ui";
char fmpz_clog_name[] = "fmpz_clog";
char fmpz_clog_ui_name[] = "fmpz_clog_ui";
char fmpz_cmp2abs_name[] = "fmpz_cmp2abs";
char fmpz_cmpabs_name[] = "fmpz_cmpabs";
char fmpz_cmp_name[] = "fmpz_cmp";
char fmpz_cmp_si_name[] = "fmpz_cmp_si";
char fmpz_cmp_ui_name[] = "fmpz_cmp_ui";
char fmpz_comb_init_clear_name[] = "fmpz_comb_init_clear";
char fmpz_combit_name[] = "fmpz_combit";
char fmpz_complement_name[] = "fmpz_complement";
char fmpz_CRT_name[] = "fmpz_CRT";
char fmpz_CRT_ui_name[] = "fmpz_CRT_ui";
char fmpz_divexact2_uiui_name[] = "fmpz_divexact2_uiui";
char fmpz_divexact_name[] = "fmpz_divexact";
char fmpz_divexact_si_name[] = "fmpz_divexact_si";
char fmpz_divexact_ui_name[] = "fmpz_divexact_ui";
char fmpz_divides_name[] = "fmpz_divides";
char fmpz_divides_mod_list_name[] = "fmpz_divides_mod_list";
char fmpz_divisible_name[] = "fmpz_divisible";
char fmpz_divisible_si_name[] = "fmpz_divisible_si";
char fmpz_divisor_in_residue_class_lenstra_name[] = "fmpz_divisor_in_residue_class_lenstra";
char fmpz_divisor_sigma_name[] = "fmpz_divisor_sigma";
char fmpz_div_newton_name[] = "fmpz_div_newton";
char fmpz_dlog_name[] = "fmpz_dlog";
char fmpz_equal_name[] = "fmpz_equal";
char fmpz_equal_si_name[] = "fmpz_equal_si";
char fmpz_equal_ui_name[] = "fmpz_equal_ui";
char fmpz_euler_phi_name[] = "fmpz_euler_phi";
char fmpz_fac_ui_name[] = "fmpz_fac_ui";
char fmpz_fdiv_q_2exp_name[] = "fmpz_fdiv_q_2exp";
char fmpz_fdiv_q_name[] = "fmpz_fdiv_q";
char fmpz_fdiv_qr_name[] = "fmpz_fdiv_qr";
char fmpz_fdiv_qr_preinvn_name[] = "fmpz_fdiv_qr_preinvn";
char fmpz_fdiv_q_si_name[] = "fmpz_fdiv_q_si";
char fmpz_fdiv_q_ui_name[] = "fmpz_fdiv_q_ui";
char fmpz_fdiv_r_2exp_name[] = "fmpz_fdiv_r_2exp";
char fmpz_fdiv_r_name[] = "fmpz_fdiv_r";
char fmpz_fdiv_ui_name[] = "fmpz_fdiv_ui";
char fmpz_fib_ui_name[] = "fmpz_fib_ui";
char fmpz_fits_si_name[] = "fmpz_fits_si";
char fmpz_flog_name[] = "fmpz_flog";
char fmpz_flog_ui_name[] = "fmpz_flog_ui";
char fmpz_fmma_name[] = "fmpz_fmma";
char fmpz_fmms_name[] = "fmpz_fmms";
char fmpz_fmpz_name[] = "fmpz_fmpz";
char fmpz_cleanup_name[] = "fmpz_cleanup";
char fmpz_stress_name[] = "fmpz_stress";
char fmpz_gcd3_name[] = "fmpz_gcd3";
char fmpz_gcd_name[] = "fmpz_gcd";
char fmpz_gcdinv_name[] = "fmpz_gcdinv";
char fmpz_gcd_ui_name[] = "fmpz_gcd_ui";
char fmpz_get_d_2exp_name[] = "fmpz_get_d_2exp";
char fmpz_get_d_name[] = "fmpz_get_d";
char fmpz_get_mpf_name[] = "fmpz_get_mpf";
char fmpz_get_mpfr_name[] = "fmpz_get_mpfr";
char fmpz_get_mpn_name[] = "fmpz_get_mpn";
char fmpz_get_mpz_name[] = "fmpz_get_mpz";
char fmpz_get_nmod_name[] = "fmpz_get_nmod";
char fmpz_get_set_ui_array_name[] = "fmpz_get_set_ui_array";
char fmpz_get_si_name[] = "fmpz_get_si";
char fmpz_get_str_name[] = "fmpz_get_str";
char fmpz_get_ui_name[] = "fmpz_get_ui";
char fmpz_init2_name[] = "fmpz_init2";
char fmpz_init_set_name[] = "fmpz_init_set";
char fmpz_init_set_readonly_name[] = "fmpz_init_set_readonly";
char fmpz_init_set_ui_name[] = "fmpz_init_set_ui";
char fmpz_invmod_name[] = "fmpz_invmod";
char fmpz_is_even_name[] = "fmpz_is_even";
char fmpz_is_perfect_power_name[] = "fmpz_is_perfect_power";
char fmpz_is_prime_name[] = "fmpz_is_prime";
char fmpz_is_prime_morrison_name[] = "fmpz_is_prime_morrison";
char fmpz_is_prime_pocklington_name[] = "fmpz_is_prime_pocklington";
char fmpz_is_prime_pseudosquare_name[] = "fmpz_is_prime_pseudosquare";
char fmpz_is_probabprime_BPSW_name[] = "fmpz_is_probabprime_BPSW";
char fmpz_is_probabprime_lucas_name[] = "fmpz_is_probabprime_lucas";
char fmpz_is_square_name[] = "fmpz_is_square";
char fmpz_is_strong_probabprime_name[] = "fmpz_is_strong_probabprime";
char fmpz_jacobi_name[] = "fmpz_jacobi";
char fmpz_kronecker_name[] = "fmpz_kronecker";
char fmpz_lcm_name[] = "fmpz_lcm";
char fmpz_mod_name[] = "fmpz_mod";
char fmpz_mod_ui_name[] = "fmpz_mod_ui";
char fmpz_moebius_mu_name[] = "fmpz_moebius_mu";
char fmpz_mpz_init_set_readonly_name[] = "fmpz_mpz_init_set_readonly";
char fmpz_mul_2exp_name[] = "fmpz_mul_2exp";
char fmpz_mul2_uiui_name[] = "fmpz_mul2_uiui";
char fmpz_mul_name[] = "fmpz_mul";
char fmpz_mul_si_name[] = "fmpz_mul_si";
char fmpz_mul_si_tdiv_q_2exp_name[] = "fmpz_mul_si_tdiv_q_2exp";
char fmpz_mul_tdiv_q_2exp_name[] = "fmpz_mul_tdiv_q_2exp";
char fmpz_multi_CRT_multi_mod_name[] = "fmpz_multi_CRT_multi_mod";
char fmpz_multi_CRT_ui_name[] = "fmpz_multi_CRT_ui";
char fmpz_mul_ui_name[] = "fmpz_mul_ui";
char fmpz_ndiv_qr_name[] = "fmpz_ndiv_qr";
char fmpz_neg_name[] = "fmpz_neg";
char fmpz_neg_ui_name[] = "fmpz_neg_ui";
char fmpz_neg_uiui_name[] = "fmpz_neg_uiui";
char fmpz_nextprime_name[] = "fmpz_nextprime";
char fmpz_or_name[] = "fmpz_or";
char fmpz_out_inp_raw_name[] = "fmpz_out_inp_raw";
char fmpz_popcnt_name[] = "fmpz_popcnt";
char fmpz_powm_name[] = "fmpz_powm";
char fmpz_powm_ui_name[] = "fmpz_powm_ui";
char fmpz_pow_ui_name[] = "fmpz_pow_ui";
char fmpz_primorial_name[] = "fmpz_primorial";
char fmpz_print_read_name[] = "fmpz_print_read";
char fmpz_randprime_name[] = "fmpz_randprime";
char fmpz_remove_name[] = "fmpz_remove";
char fmpz_rfac_ui_name[] = "fmpz_rfac_ui";
char fmpz_rfac_uiui_name[] = "fmpz_rfac_uiui";
char fmpz_root_name[] = "fmpz_root";
char fmpz_setbit_name[] = "fmpz_setbit";
char fmpz_set_name[] = "fmpz_set";
char fmpz_set_d_2exp_name[] = "fmpz_set_d_2exp";
char fmpz_set_signed_ui_array_name[] = "fmpz_set_signed_ui_array";
char fmpz_set_signed_uiui_name[] = "fmpz_set_signed_uiui";
char fmpz_set_signed_uiuiui_name[] = "fmpz_set_signed_uiuiui";
char fmpz_set_str_name[] = "fmpz_set_str";
char fmpz_set_ui_smod_name[] = "fmpz_set_ui_smod";
char fmpz_set_uiui_name[] = "fmpz_set_uiui";
char fmpz_sgn_name[] = "fmpz_sgn";
char fmpz_size_name[] = "fmpz_size";
char fmpz_sizeinbase_name[] = "fmpz_sizeinbase";
char fmpz_smod_name[] = "fmpz_smod";
char fmpz_sqrt_name[] = "fmpz_sqrt";
char fmpz_sqrtmod_name[] = "fmpz_sqrtmod";
char fmpz_sqrtrem_name[] = "fmpz_sqrtrem";
char fmpz_sub_name[] = "fmpz_sub";
char fmpz_submul_name[] = "fmpz_submul";
char fmpz_submul_si_name[] = "fmpz_submul_si";
char fmpz_submul_ui_name[] = "fmpz_submul_ui";
char fmpz_sub_ui_name[] = "fmpz_sub_ui";
char fmpz_swap_name[] = "fmpz_swap";
char fmpz_tdiv_q_2exp_name[] = "fmpz_tdiv_q_2exp";
char fmpz_tdiv_q_name[] = "fmpz_tdiv_q";
char fmpz_tdiv_qr_name[] = "fmpz_tdiv_qr";
char fmpz_tdiv_q_si_name[] = "fmpz_tdiv_q_si";
char fmpz_tdiv_q_ui_name[] = "fmpz_tdiv_q_ui";
char fmpz_tdiv_r_2exp_name[] = "fmpz_tdiv_r_2exp";
char fmpz_tdiv_ui_name[] = "fmpz_tdiv_ui";
char fmpz_tstbit_name[] = "fmpz_tstbit";
char fmpz_val2_name[] = "fmpz_val2";
char fmpz_xgcd_name[] = "fmpz_xgcd";
char fmpz_xgcd_canonical_bezout_name[] = "fmpz_xgcd_canonical_bezout";
char fmpz_xgcd_partial_name[] = "fmpz_xgcd_partial";
char fmpz_xor_name[] = "fmpz_xor";

char * test_names[] =
{
    fmpz_abs_name,
    fmpz_abs_fits_ui_name,
    fmpz_abs_lbound_ui_2exp_name,
    fmpz_abs_ubound_ui_2exp_name,
    fmpz_add_name,
    fmpz_addmul_name,
    fmpz_addmul_si_name,
    fmpz_addmul_ui_name,
    fmpz_add_ui_name,
    fmpz_and_name,
    fmpz_bin_uiui_name,
    fmpz_bit_pack_name,
    fmpz_bits_name,
    fmpz_cdiv_q_2exp_name,
    fmpz_cdiv_q_name,
    fmpz_cdiv_qr_name,
    fmpz_cdiv_q_si_name,
    fmpz_cdiv_q_ui_name,
    fmpz_cdiv_r_2exp_name,
    fmpz_cdiv_ui_name,
    fmpz_clog_name,
    fmpz_clog_ui_name,
    fmpz_cmp2abs_name,
    fmpz_cmpabs_name,
    fmpz_cmp_name,
    fmpz_cmp_si_name,
    fmpz_cmp_ui_name,
    fmpz_comb_init_clear_name,
    fmpz_combit_name,
    fmpz_complement_name,
    fmpz_CRT_name,
    fmpz_CRT_ui_name,
    fmpz_divexact2_uiui_name,
    fmpz_divexact_name,
    fmpz_divexact_si_name,
    fmpz_divexact_ui_name,
    fmpz_divides_name,
    fmpz_divides_mod_list_name,
    fmpz_divisible_name,
    fmpz_divisible_si_name,
    fmpz_divisor_in_residue_class_lenstra_name,
    fmpz_divisor_sigma_name,
    fmpz_div_newton_name,
    fmpz_dlog_name,
    fmpz_equal_name,
    fmpz_equal_si_name,
    fmpz_equal_ui_name,
    fmpz_euler_phi_name,
    fmpz_fac_ui_name,
    fmpz_fdiv_q_2exp_name,
    fmpz_fdiv_q_name,
    fmpz_fdiv_qr_name,
    fmpz_fdiv_qr_preinvn_name,
    fmpz_fdiv_q_si_name,
    fmpz_fdiv_q_ui_name,
    fmpz_fdiv_r_2exp_name,
    fmpz_fdiv_r_name,
    fmpz_fdiv_ui_name,
    fmpz_fib_ui_name,
    fmpz_fits_si_name,
    fmpz_flog_name,
    fmpz_flog_ui_name,
    fmpz_fmma_name,
    fmpz_fmms_name,
    fmpz_fmpz_name,
    fmpz_cleanup_name,
    fmpz_stress_name,
    fmpz_gcd3_name,
    fmpz_gcd_name,
    fmpz_gcdinv_name,
    fmpz_gcd_ui_name,
    fmpz_get_d_2exp_name,
    fmpz_get_d_name,
    fmpz_get_mpf_name,
    fmpz_get_mpfr_name,
    fmpz_get_mpn_name,
    fmpz_get_mpz_name,
    fmpz_get_nmod_name,
    fmpz_get_set_ui_array_name,
    fmpz_get_si_name,
    fmpz_get_str_name,
    fmpz_get_ui_name,
    fmpz_init2_name,
    fmpz_init_set_name,
    fmpz_init_set_readonly_name,
    fmpz_init_set_ui_name,
    fmpz_invmod_name,
    fmpz_is_even_name,
    fmpz_is_perfect_power_name,
    fmpz_is_prime_name,
    fmpz_is_prime_morrison_name,
    fmpz_is_prime_pocklington_name,
    fmpz_is_prime_pseudosquare_name,
    fmpz_is_probabprime_BPSW_name,
    fmpz_is_probabprime_lucas_name,
    fmpz_is_square_name,
    fmpz_is_strong_probabprime_name,
    fmpz_jacobi_name,
    fmpz_kronecker_name,
    fmpz_lcm_name,
    fmpz_mod_name,
    fmpz_mod_ui_name,
    fmpz_moebius_mu_name,
    fmpz_mpz_init_set_readonly_name,
    fmpz_mul_2exp_name,
    fmpz_mul2_uiui_name,
    fmpz_mul_name,
    fmpz_mul_si_name,
    fmpz_mul_si_tdiv_q_2exp_name,
    fmpz_mul_tdiv_q_2exp_name,
    fmpz_multi_CRT_multi_mod_name,
    fmpz_multi_CRT_ui_name,
    fmpz_mul_ui_name,
    fmpz_ndiv_qr_name,
    fmpz_neg_name,
    fmpz_neg_ui_name,
    fmpz_neg_uiui_name,
    fmpz_nextprime_name,
    fmpz_or_name,
    fmpz_out_inp_raw_name,
    fmpz_popcnt_name,
    fmpz_powm_name,
    fmpz_powm_ui_name,
    fmpz_pow_ui_name,
    fmpz_primorial_name,
    fmpz_print_read_name,
    fmpz_randprime_name,
    fmpz_remove_name,
    fmpz_rfac_ui_name,
    fmpz_rfac_uiui_name,
    fmpz_root_name,
    fmpz_setbit_name,
    fmpz_set_name,
    fmpz_set_d_2exp_name,
    fmpz_set_signed_ui_array_name,
    fmpz_set_signed_uiui_name,
    fmpz_set_signed_uiuiui_name,
    fmpz_set_str_name,
    fmpz_set_ui_smod_name,
    fmpz_set_uiui_name,
    fmpz_sgn_name,
    fmpz_size_name,
    fmpz_sizeinbase_name,
    fmpz_smod_name,
    fmpz_sqrt_name,
    fmpz_sqrtmod_name,
    fmpz_sqrtrem_name,
    fmpz_sub_name,
    fmpz_submul_name,
    fmpz_submul_si_name,
    fmpz_submul_ui_name,
    fmpz_sub_ui_name,
    fmpz_swap_name,
    fmpz_tdiv_q_2exp_name,
    fmpz_tdiv_q_name,
    fmpz_tdiv_qr_name,
    fmpz_tdiv_q_si_name,
    fmpz_tdiv_q_ui_name,
    fmpz_tdiv_r_2exp_name,
    fmpz_tdiv_ui_name,
    fmpz_tstbit_name,
    fmpz_val2_name,
    fmpz_xgcd_name,
    fmpz_xgcd_canonical_bezout_name,
    fmpz_xgcd_partial_name,
    fmpz_xor_name
};

/* main function *************************************************************/

#define NUMBER_OF_TESTS (sizeof(test_functions) / sizeof(int (*)(void)))

int
main(int argc, char ** argv)
{
    int ix, jx;

    if (argc < 2)
    {
        for (ix = 0; ix < NUMBER_OF_TESTS; ix++)
            if (test_functions[ix]())
                flint_abort();
    }
    else
    {
        for (ix = 1; ix < argc; ix++)
        {
            for (jx = 0; jx < NUMBER_OF_TESTS; jx++)
            {
                /* If argument equals to test name, run it */
                if (strcmp(argv[ix], test_names[jx]) == 0)
                {
                    if (test_functions[jx]())
                        flint_abort();
                    break;
                }
            }

            if (jx == NUMBER_OF_TESTS)
            {
                fprintf(stderr, "Error: Could not find test function for %s\n", argv[ix]);
                flint_abort();
            }
        }
    }

    return 0;
}
