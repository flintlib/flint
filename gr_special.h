/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef GR_SPECIAL_H
#define GR_SPECIAL_H

#ifdef GR_SPECIAL_INLINES_C
#define GR_SPECIAL_INLINE FLINT_DLL
#else
#define GR_SPECIAL_INLINE static __inline__
#endif

#include "gr.h"

#ifdef __cplusplus
 extern "C" {
#endif

#define GR_SPECIAL_DEF GR_SPECIAL_INLINE WARN_UNUSED_RESULT

GR_SPECIAL_DEF int gr_pi(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, PI)(res, ctx); }
GR_SPECIAL_DEF int gr_euler(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, EULER)(res, ctx); }
GR_SPECIAL_DEF int gr_catalan(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, CATALAN)(res, ctx); }
GR_SPECIAL_DEF int gr_khinchin(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, KHINCHIN)(res, ctx); }
GR_SPECIAL_DEF int gr_glaisher(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, GLAISHER)(res, ctx); }

GR_SPECIAL_DEF int gr_exp(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, EXP)(res, x, ctx); }
GR_SPECIAL_DEF int gr_expm1(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, EXPM1)(res, x, ctx); }
GR_SPECIAL_DEF int gr_exp2(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, EXP2)(res, x, ctx); }
GR_SPECIAL_DEF int gr_exp10(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, EXP10)(res, x, ctx); }
GR_SPECIAL_DEF int gr_exp_pi_i(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, EXP_PI_I)(res, x, ctx); }
GR_SPECIAL_DEF int gr_log(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, LOG)(res, x, ctx); }
GR_SPECIAL_DEF int gr_log1p(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, LOG1P)(res, x, ctx); }
GR_SPECIAL_DEF int gr_log_pi_i(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, LOG_PI_I)(res, x, ctx); }
GR_SPECIAL_DEF int gr_log2(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, LOG2)(res, x, ctx); }
GR_SPECIAL_DEF int gr_log10(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, LOG10)(res, x, ctx); }

GR_SPECIAL_DEF int gr_sin(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SIN)(res, x, ctx); }
GR_SPECIAL_DEF int gr_cos(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, COS)(res, x, ctx); }
GR_SPECIAL_DEF int gr_sin_cos(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_ctx_t ctx) { return GR_BINARY_UNARY_OP(ctx, SIN_COS)(res1, res2, x, ctx); }
GR_SPECIAL_DEF int gr_tan(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, TAN)(res, x, ctx); }
GR_SPECIAL_DEF int gr_cot(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, COT)(res, x, ctx); }
GR_SPECIAL_DEF int gr_sec(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SEC)(res, x, ctx); }
GR_SPECIAL_DEF int gr_csc(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, CSC)(res, x, ctx); }

GR_SPECIAL_DEF int gr_sin_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SIN_PI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_cos_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, COS_PI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_sin_cos_pi(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_ctx_t ctx) { return GR_BINARY_UNARY_OP(ctx, SIN_COS_PI)(res1, res2, x, ctx); }
GR_SPECIAL_DEF int gr_tan_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, TAN_PI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_cot_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, COT_PI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_sec_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SEC_PI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_csc_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, CSC_PI)(res, x, ctx); }

GR_SPECIAL_DEF int gr_sinc(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SINC)(res, x, ctx); }
GR_SPECIAL_DEF int gr_sinc_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SINC_PI)(res, x, ctx); }

GR_SPECIAL_DEF int gr_sinh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SINH)(res, x, ctx); }
GR_SPECIAL_DEF int gr_cosh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, COSH)(res, x, ctx); }
GR_SPECIAL_DEF int gr_sinh_cosh(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_ctx_t ctx) { return GR_BINARY_UNARY_OP(ctx, SINH_COSH)(res1, res2, x, ctx); }
GR_SPECIAL_DEF int gr_tanh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, TANH)(res, x, ctx); }
GR_SPECIAL_DEF int gr_coth(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, COTH)(res, x, ctx); }
GR_SPECIAL_DEF int gr_sech(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SECH)(res, x, ctx); }
GR_SPECIAL_DEF int gr_csch(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, CSCH)(res, x, ctx); }

GR_SPECIAL_DEF int gr_asin(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ASIN)(res, x, ctx); }
GR_SPECIAL_DEF int gr_acos(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ACOS)(res, x, ctx); }
GR_SPECIAL_DEF int gr_atan(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ATAN)(res, x, ctx); }
GR_SPECIAL_DEF int gr_atan2(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, ATAN2)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_acot(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ACOT)(res, x, ctx); }
GR_SPECIAL_DEF int gr_asec(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ASEC)(res, x, ctx); }
GR_SPECIAL_DEF int gr_acsc(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ACSC)(res, x, ctx); }

GR_SPECIAL_DEF int gr_asin_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ASIN_PI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_acos_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ACOS_PI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_atan_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ATAN_PI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_acot_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ACOT_PI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_asec_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ASEC_PI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_acsc_pi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ACSC_PI)(res, x, ctx); }

GR_SPECIAL_DEF int gr_asinh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ASINH)(res, x, ctx); }
GR_SPECIAL_DEF int gr_acosh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ACOSH)(res, x, ctx); }
GR_SPECIAL_DEF int gr_atanh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ATANH)(res, x, ctx); }
GR_SPECIAL_DEF int gr_acoth(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ACOTH)(res, x, ctx); }
GR_SPECIAL_DEF int gr_asech(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ASECH)(res, x, ctx); }
GR_SPECIAL_DEF int gr_acsch(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ACSCH)(res, x, ctx); }

GR_SPECIAL_DEF int gr_fac(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, FAC)(res, x, ctx); }
GR_SPECIAL_DEF int gr_fac_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_UI(ctx, FAC_UI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_fac_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx) { return GR_UNARY_OP_FMPZ(ctx, FAC_FMPZ)(res, x, ctx); }
GR_SPECIAL_DEF int gr_fac_vec(gr_ptr res, slong len, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, FAC_VEC)(res, len, ctx); }

GR_SPECIAL_DEF int gr_rfac_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_UI(ctx, RFAC_UI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_rfac_vec(gr_ptr res, slong len, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, RFAC_VEC)(res, len, ctx); }

GR_SPECIAL_DEF int gr_bin(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, BIN)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_bin_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, BIN_UI)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_bin_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI_UI(ctx, BIN_UI_UI)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_bin_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx) { return GR_BINARY_OP_UI_SI(ctx, BIN_UI_VEC)(res, x, len, ctx); }

GR_SPECIAL_DEF int gr_rising(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, RISING)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_rising_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, RISING_UI)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_falling(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, FALLING)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_falling_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, FALLING_UI)(res, x, y, ctx); }

GR_SPECIAL_DEF int gr_gamma(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, GAMMA)(res, x, ctx); }
GR_SPECIAL_DEF int gr_gamma_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx) { return GR_UNARY_OP_FMPZ(ctx, GAMMA_FMPZ)(res, x, ctx); }
GR_SPECIAL_DEF int gr_gamma_fmpq(gr_ptr res, const fmpq_t x, gr_ctx_t ctx) { return GR_UNARY_OP_FMPQ(ctx, GAMMA_FMPQ)(res, x, ctx); }
GR_SPECIAL_DEF int gr_rgamma(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, RGAMMA)(res, x, ctx); }
GR_SPECIAL_DEF int gr_lgamma(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, LGAMMA)(res, x, ctx); }
GR_SPECIAL_DEF int gr_digamma(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, DIGAMMA)(res, x, ctx); }
GR_SPECIAL_DEF int gr_beta(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, BETA)(res, x, y, ctx); }

GR_SPECIAL_DEF int gr_bernoulli_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_UI(ctx, BERNOULLI_UI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_bernoulli_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx) { return GR_UNARY_OP_FMPZ(ctx, BERNOULLI_FMPZ)(res, x, ctx); }
GR_SPECIAL_DEF int gr_bernoulli_vec(gr_ptr res, slong len, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, BERNOULLI_VEC)(res, len, ctx); }

GR_SPECIAL_DEF int gr_harmonic_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_UI(ctx, HARMONIC_UI)(res, x, ctx); }

GR_SPECIAL_DEF int gr_fib_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_UI(ctx, FIB_UI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_fib_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx) { return GR_UNARY_OP_FMPZ(ctx, FIB_FMPZ)(res, x, ctx); }
GR_SPECIAL_DEF int gr_fib_vec(gr_ptr res, slong len, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, FIB_VEC)(res, len, ctx); }

GR_SPECIAL_DEF int gr_eulernum_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_UI(ctx, EULERNUM_UI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_eulernum_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx) { return GR_UNARY_OP_FMPZ(ctx, EULERNUM_FMPZ)(res, x, ctx); }
GR_SPECIAL_DEF int gr_eulernum_vec(gr_ptr res, slong len, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, EULERNUM_VEC)(res, len, ctx); }

GR_SPECIAL_DEF int gr_bellnum_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_UI(ctx, BELLNUM_UI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_bellnum_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx) { return GR_UNARY_OP_FMPZ(ctx, BELLNUM_FMPZ)(res, x, ctx); }
GR_SPECIAL_DEF int gr_bellnum_vec(gr_ptr res, slong len, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, BELLNUM_VEC)(res, len, ctx); }

GR_SPECIAL_DEF int gr_doublefac_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_UI(ctx, DOUBLEFAC_UI)(res, x, ctx); }

GR_SPECIAL_DEF int gr_stirling_s1u_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI_UI(ctx, STIRLING_S1U_UIUI)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_stirling_s1_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI_UI(ctx, STIRLING_S1_UIUI)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_stirling_s2_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI_UI(ctx, STIRLING_S2_UIUI)(res, x, y, ctx); }

GR_SPECIAL_DEF int gr_stirling_s1u_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx) { return GR_BINARY_OP_UI_SI(ctx, STIRLING_S1U_UI_VEC)(res, x, len, ctx); }
GR_SPECIAL_DEF int gr_stirling_s1_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx) { return GR_BINARY_OP_UI_SI(ctx, STIRLING_S1_UI_VEC)(res, x, len, ctx); }
GR_SPECIAL_DEF int gr_stirling_s2_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx) { return GR_BINARY_OP_UI_SI(ctx, STIRLING_S2_UI_VEC)(res, x, len, ctx); }

GR_SPECIAL_DEF int gr_partitions_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_UI(ctx, PARTITIONS_UI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_partitions_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx) { return GR_UNARY_OP_FMPZ(ctx, PARTITIONS_FMPZ)(res, x, ctx); }
GR_SPECIAL_DEF int gr_partitions_vec(gr_ptr res, slong len, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, PARTITIONS_VEC)(res, len, ctx); }

GR_SPECIAL_DEF int gr_erf(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ERF)(res, x, ctx); }
GR_SPECIAL_DEF int gr_erfc(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ERFC)(res, x, ctx); }
GR_SPECIAL_DEF int gr_erfcx(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ERFCX)(res, x, ctx); }
GR_SPECIAL_DEF int gr_erfi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ERFI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_erfinv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ERFINV)(res, x, ctx); }
GR_SPECIAL_DEF int gr_erfcinv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ERFCINV)(res, x, ctx); }

/* todo: normalized flags */
GR_SPECIAL_DEF int gr_fresnel_s(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, FRESNEL_S)(res, x, ctx); }
GR_SPECIAL_DEF int gr_fresnel_c(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, FRESNEL_C)(res, x, ctx); }
GR_SPECIAL_DEF int gr_fresnel(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_ctx_t ctx) { return GR_BINARY_UNARY_OP(ctx, FRESNEL)(res1, res2, x, ctx); }
GR_SPECIAL_DEF int gr_gamma_upper(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, GAMMA_UPPER)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_gamma_lower(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, GAMMA_LOWER)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_beta_lower(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ctx_t ctx) { return GR_TERNARY_OP(ctx, BETA_LOWER)(res, x, y, z, ctx); }
GR_SPECIAL_DEF int gr_exp_integral(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, EXP_INTEGRAL)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_exp_integral_ei(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, EXP_INTEGRAL_EI)(res, x, ctx); }

GR_SPECIAL_DEF int gr_sin_integral(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SIN_INTEGRAL)(res, x, ctx); }
GR_SPECIAL_DEF int gr_cos_integral(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, COS_INTEGRAL)(res, x, ctx); }
GR_SPECIAL_DEF int gr_sinh_integral(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SINH_INTEGRAL)(res, x, ctx); }
GR_SPECIAL_DEF int gr_cosh_integral(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, COSH_INTEGRAL)(res, x, ctx); }
GR_SPECIAL_DEF int gr_log_integral(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, LOG_INTEGRAL)(res, x, ctx); }
GR_SPECIAL_DEF int gr_dilog(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, DILOG)(res, x, ctx); }

GR_SPECIAL_DEF int gr_zeta(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ZETA)(res, x, ctx); }
GR_SPECIAL_DEF int gr_zeta_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_UI(ctx, ZETA_UI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_hurwitz_zeta(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, HURWITZ_ZETA)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_lerch_phi(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ctx_t ctx) { return GR_TERNARY_OP(ctx, LERCH_PHI)(res, x, y, z, ctx); }
GR_SPECIAL_DEF int gr_stieltjes(gr_ptr res, const fmpz_t x, gr_srcptr y, gr_ctx_t ctx) { return GR_FMPZ_BINARY_OP(ctx, STIELTJES)(res, x, y, ctx); }



/* todo

    GR_METHOD_CHEBYSHEV_T_UI,
    GR_METHOD_CHEBYSHEV_T_FMPZ,
    GR_METHOD_CHEBYSHEV_T,
    GR_METHOD_CHEBYSHEV_U_UI,
    GR_METHOD_CHEBYSHEV_U_FMPZ,
    GR_METHOD_CHEBYSHEV_U,
    GR_METHOD_JACOBI_P,
    GR_METHOD_GEGENBAUER_C,
    GR_METHOD_LAGUERRE_L,
    GR_METHOD_HERMITE_H,
    GR_METHOD_LEGENDRE_P,
    GR_METHOD_LEGENDRE_Q,
    GR_METHOD_GAUSS_LEGENDRE_NODE,
    GR_METHOD_SPHERICAL_Y_SI,

    GR_METHOD_LAMBERTW,
    GR_METHOD_LAMBERTW_FMPZ,

    GR_METHOD_BESSEL_J,
    GR_METHOD_BESSEL_Y,
    GR_METHOD_BESSEL_J_Y,
    GR_METHOD_BESSEL_I,
    GR_METHOD_BESSEL_I_SCALED,
    GR_METHOD_BESSEL_K,
    GR_METHOD_BESSEL_K_SCALED,

    GR_METHOD_AIRY,
    GR_METHOD_AIRY_AI,
    GR_METHOD_AIRY_BI,
    GR_METHOD_AIRY_AI_PRIME,
    GR_METHOD_AIRY_BI_PRIME,
    GR_METHOD_AIRY_AI_ZERO_FMPZ,
    GR_METHOD_AIRY_BI_ZERO_FMPZ,
    GR_METHOD_AIRY_AI_PRIME_ZERO_FMPZ,
    GR_METHOD_AIRY_BI_PRIME_ZERO_FMPZ,

    GR_METHOD_COULOMB_F,
    GR_METHOD_COULOMB_G,
    GR_METHOD_COULOMB,

    GR_METHOD_DIRICHLET_CHI_UI,
    GR_METHOD_DIRICHLET_L,

    GR_METHOD_BERNPOLY_UI,
    GR_METHOD_EULERPOLY_UI,
    GR_METHOD_POLYLOG,
    GR_METHOD_POLYLOG_SI,
    GR_METHOD_BARNESG,
    GR_METHOD_LOG_BARNESG,
    GR_METHOD_POLYGAMMA,

    GR_METHOD_HYPGEOM_0F1,
    GR_METHOD_HYPGEOM_1F1,
    GR_METHOD_HYPGEOM_2F0,
    GR_METHOD_HYPGEOM_2F1,
    GR_METHOD_HYPGEOM_U,
    GR_METHOD_HYPGEOM_PFQ,

    GR_METHOD_AGM,
    GR_METHOD_AGM1,

*/

/* todo: elliptic, modular */


/* generic implementations */

int gr_generic_fib2_fmpz(gr_ptr v, gr_ptr u, const fmpz_t n, gr_ctx_t ctx);
int gr_generic_fib_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx);
int gr_generic_fib_ui(gr_ptr res, ulong n, gr_ctx_t ctx);
int gr_generic_fib_vec(gr_ptr res, slong len, gr_ctx_t ctx);


#ifdef __cplusplus
}
#endif

#endif
