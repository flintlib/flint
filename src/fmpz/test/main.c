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
#include <stdlib.h>
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
#include "t-and.c"
#include "t-aors_ui.c"
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

test_struct tests[] =
{
    TEST_FUNCTION(fmpz_abs),
    TEST_FUNCTION(fmpz_abs_fits_ui),
    TEST_FUNCTION(fmpz_abs_lbound_ui_2exp),
    TEST_FUNCTION(fmpz_abs_ubound_ui_2exp),
    TEST_FUNCTION(fmpz_add),
    TEST_FUNCTION(fmpz_addmul),
    TEST_FUNCTION(fmpz_addmul_si),
    TEST_FUNCTION(fmpz_addmul_ui),
    TEST_FUNCTION(fmpz_and),
    TEST_FUNCTION(fmpz_aors_ui),
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

/* main function *************************************************************/

TEST_MAIN(tests)
