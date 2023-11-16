/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_SPECIAL_H
#define GR_SPECIAL_H

#ifdef GR_SPECIAL_INLINES_C
#define GR_SPECIAL_INLINE
#else
#define GR_SPECIAL_INLINE static inline
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

GR_SPECIAL_DEF int gr_lambertw(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, LAMBERTW)(res, x, ctx); }
GR_SPECIAL_DEF int gr_lambertw_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t k, gr_ctx_t ctx) { return GR_BINARY_OP_FMPZ(ctx, LAMBERTW_FMPZ)(res, x, k, ctx); }

GR_SPECIAL_DEF int gr_fac(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, FAC)(res, x, ctx); }
GR_SPECIAL_DEF int gr_fac_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_UI(ctx, FAC_UI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_fac_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx) { return GR_UNARY_OP_FMPZ(ctx, FAC_FMPZ)(res, x, ctx); }
GR_SPECIAL_DEF int gr_fac_vec(gr_ptr res, slong len, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, FAC_VEC)(res, len, ctx); }

GR_SPECIAL_DEF int gr_rfac(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, RFAC)(res, x, ctx); }
GR_SPECIAL_DEF int gr_rfac_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_UI(ctx, RFAC_UI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_rfac_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx) { return GR_UNARY_OP_FMPZ(ctx, RFAC_FMPZ)(res, x, ctx); }
GR_SPECIAL_DEF int gr_rfac_vec(gr_ptr res, slong len, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, RFAC_VEC)(res, len, ctx); }

GR_SPECIAL_DEF int gr_bin(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, BIN)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_bin_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, BIN_UI)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_bin_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI_UI(ctx, BIN_UIUI)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_bin_vec(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, BIN_VEC)(res, x, len, ctx); }
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

GR_SPECIAL_DEF int gr_barnes_g(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, BARNES_G)(res, x, ctx); }
GR_SPECIAL_DEF int gr_log_barnes_g(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, LOG_BARNES_G)(res, x, ctx); }

GR_SPECIAL_DEF int gr_doublefac(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, DOUBLEFAC)(res, x, ctx); }
GR_SPECIAL_DEF int gr_doublefac_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_UI(ctx, DOUBLEFAC_UI)(res, x, ctx); }

GR_SPECIAL_DEF int gr_harmonic(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, HARMONIC)(res, x, ctx); }
GR_SPECIAL_DEF int gr_harmonic_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_UI(ctx, HARMONIC_UI)(res, x, ctx); }

GR_SPECIAL_DEF int gr_bernoulli_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_UI(ctx, BERNOULLI_UI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_bernoulli_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx) { return GR_UNARY_OP_FMPZ(ctx, BERNOULLI_FMPZ)(res, x, ctx); }
GR_SPECIAL_DEF int gr_bernoulli_vec(gr_ptr res, slong len, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, BERNOULLI_VEC)(res, len, ctx); }

GR_SPECIAL_DEF int gr_fib_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_UI(ctx, FIB_UI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_fib_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx) { return GR_UNARY_OP_FMPZ(ctx, FIB_FMPZ)(res, x, ctx); }
GR_SPECIAL_DEF int gr_fib_vec(gr_ptr res, slong len, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, FIB_VEC)(res, len, ctx); }

GR_SPECIAL_DEF int gr_eulernum_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_UI(ctx, EULERNUM_UI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_eulernum_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx) { return GR_UNARY_OP_FMPZ(ctx, EULERNUM_FMPZ)(res, x, ctx); }
GR_SPECIAL_DEF int gr_eulernum_vec(gr_ptr res, slong len, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, EULERNUM_VEC)(res, len, ctx); }

GR_SPECIAL_DEF int gr_bernpoly_ui(gr_ptr res, ulong n, gr_srcptr x, gr_ctx_t ctx) { return GR_UI_BINARY_OP(ctx, BERNPOLY_UI)(res, n, x, ctx); }
GR_SPECIAL_DEF int gr_eulerpoly_ui(gr_ptr res, ulong n, gr_srcptr x, gr_ctx_t ctx) { return GR_UI_BINARY_OP(ctx, EULERPOLY_UI)(res, n, x, ctx); }

GR_SPECIAL_DEF int gr_bellnum_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_UI(ctx, BELLNUM_UI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_bellnum_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx) { return GR_UNARY_OP_FMPZ(ctx, BELLNUM_FMPZ)(res, x, ctx); }
GR_SPECIAL_DEF int gr_bellnum_vec(gr_ptr res, slong len, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, BELLNUM_VEC)(res, len, ctx); }

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
GR_SPECIAL_DEF int gr_fresnel_s(gr_ptr res, gr_srcptr x, int normalized, gr_ctx_t ctx) { return GR_UNARY_OP_WITH_FLAG(ctx, FRESNEL_S)(res, x, normalized, ctx); }
GR_SPECIAL_DEF int gr_fresnel_c(gr_ptr res, gr_srcptr x, int normalized, gr_ctx_t ctx) { return GR_UNARY_OP_WITH_FLAG(ctx, FRESNEL_C)(res, x, normalized, ctx); }
GR_SPECIAL_DEF int gr_fresnel(gr_ptr res1, gr_ptr res2, gr_srcptr x, int normalized, gr_ctx_t ctx) { return GR_BINARY_UNARY_OP_WITH_FLAG(ctx, FRESNEL)(res1, res2, x, normalized, ctx); }
GR_SPECIAL_DEF int gr_gamma_upper(gr_ptr res, gr_srcptr x, gr_srcptr y, int regularized, gr_ctx_t ctx) { return GR_BINARY_OP_WITH_FLAG(ctx, GAMMA_UPPER)(res, x, y, regularized, ctx); }
GR_SPECIAL_DEF int gr_gamma_lower(gr_ptr res, gr_srcptr x, gr_srcptr y, int regularized, gr_ctx_t ctx) { return GR_BINARY_OP_WITH_FLAG(ctx, GAMMA_LOWER)(res, x, y, regularized, ctx); }
GR_SPECIAL_DEF int gr_beta_lower(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, int regularized, gr_ctx_t ctx) { return GR_TERNARY_OP_WITH_FLAG(ctx, BETA_LOWER)(res, x, y, z, regularized, ctx); }
GR_SPECIAL_DEF int gr_exp_integral(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, EXP_INTEGRAL)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_exp_integral_ei(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, EXP_INTEGRAL_EI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_sin_integral(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SIN_INTEGRAL)(res, x, ctx); }
GR_SPECIAL_DEF int gr_cos_integral(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, COS_INTEGRAL)(res, x, ctx); }
GR_SPECIAL_DEF int gr_sinh_integral(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SINH_INTEGRAL)(res, x, ctx); }
GR_SPECIAL_DEF int gr_cosh_integral(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, COSH_INTEGRAL)(res, x, ctx); }
GR_SPECIAL_DEF int gr_log_integral(gr_ptr res, gr_srcptr x, int offset, gr_ctx_t ctx) { return GR_UNARY_OP_WITH_FLAG(ctx, LOG_INTEGRAL)(res, x, offset, ctx); }
GR_SPECIAL_DEF int gr_dilog(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, DILOG)(res, x, ctx); }

GR_SPECIAL_DEF int gr_bessel_j(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, BESSEL_J)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_bessel_y(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, BESSEL_Y)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_bessel_i(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, BESSEL_I)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_bessel_k(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, BESSEL_K)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_bessel_j_y(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_BINARY_OP(ctx, BESSEL_J_Y)(res1, res2, x, y, ctx); }
GR_SPECIAL_DEF int gr_bessel_i_scaled(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, BESSEL_I_SCALED)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_bessel_k_scaled(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, BESSEL_K_SCALED)(res, x, y, ctx); }

GR_SPECIAL_DEF int gr_airy(gr_ptr res1, gr_ptr res2, gr_ptr res3, gr_ptr res4, gr_srcptr x, gr_ctx_t ctx) { return GR_QUATERNARY_UNARY_OP(ctx, AIRY)(res1, res2, res3, res4, x, ctx); }
GR_SPECIAL_DEF int gr_airy_ai(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, AIRY_AI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_airy_bi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, AIRY_BI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_airy_ai_prime(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, AIRY_AI_PRIME)(res, x, ctx); }
GR_SPECIAL_DEF int gr_airy_bi_prime(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, AIRY_BI_PRIME)(res, x, ctx); }

GR_SPECIAL_DEF int gr_airy_ai_zero(gr_ptr res, const fmpz_t n, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, AIRY_AI_ZERO)(res, n, ctx); }
GR_SPECIAL_DEF int gr_airy_bi_zero(gr_ptr res, const fmpz_t n, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, AIRY_BI_ZERO)(res, n, ctx); }
GR_SPECIAL_DEF int gr_airy_ai_prime_zero(gr_ptr res, const fmpz_t n, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, AIRY_AI_PRIME_ZERO)(res, n, ctx); }
GR_SPECIAL_DEF int gr_airy_bi_prime_zero(gr_ptr res, const fmpz_t n, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, AIRY_BI_PRIME_ZERO)(res, n, ctx); }

GR_SPECIAL_DEF int gr_coulomb(gr_ptr res1, gr_ptr res2, gr_ptr res3, gr_ptr res4, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ctx_t ctx) { return GR_QUATERNARY_TERNARY_OP(ctx, COULOMB)(res1, res2, res3, res4, x, y, z, ctx); }
GR_SPECIAL_DEF int gr_coulomb_f(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ctx_t ctx) { return GR_TERNARY_OP(ctx, COULOMB_F)(res, x, y, z, ctx); }
GR_SPECIAL_DEF int gr_coulomb_g(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ctx_t ctx) { return GR_TERNARY_OP(ctx, COULOMB_G)(res, x, y, z, ctx); }
GR_SPECIAL_DEF int gr_coulomb_hpos(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ctx_t ctx) { return GR_TERNARY_OP(ctx, COULOMB_HPOS)(res, x, y, z, ctx); }
GR_SPECIAL_DEF int gr_coulomb_hneg(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ctx_t ctx) { return GR_TERNARY_OP(ctx, COULOMB_HNEG)(res, x, y, z, ctx); }

GR_SPECIAL_DEF int gr_chebyshev_t_fmpz(gr_ptr res, const fmpz_t n, gr_srcptr x, gr_ctx_t ctx) { return GR_FMPZ_BINARY_OP(ctx, CHEBYSHEV_T_FMPZ)(res, n, x, ctx); }
GR_SPECIAL_DEF int gr_chebyshev_t(gr_ptr res, gr_srcptr n, gr_srcptr x, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, CHEBYSHEV_T)(res, n, x, ctx); }
GR_SPECIAL_DEF int gr_chebyshev_u_fmpz(gr_ptr res, const fmpz_t n, gr_srcptr x, gr_ctx_t ctx) { return GR_FMPZ_BINARY_OP(ctx, CHEBYSHEV_U_FMPZ)(res, n, x, ctx); }
GR_SPECIAL_DEF int gr_chebyshev_u(gr_ptr res, gr_srcptr n, gr_srcptr x, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, CHEBYSHEV_U)(res, n, x, ctx); }

GR_SPECIAL_DEF int gr_jacobi_p(gr_ptr res, gr_srcptr n, gr_srcptr a, gr_srcptr b, gr_srcptr z, gr_ctx_t ctx) { return GR_QUATERNARY_OP(ctx, JACOBI_P)(res, n, a, b, z, ctx); }
GR_SPECIAL_DEF int gr_gegenbauer_c(gr_ptr res, gr_srcptr n, gr_srcptr m, gr_srcptr z, gr_ctx_t ctx) { return GR_TERNARY_OP(ctx, GEGENBAUER_C)(res, n, m, z, ctx); }
GR_SPECIAL_DEF int gr_laguerre_l(gr_ptr res, gr_srcptr n, gr_srcptr m, gr_srcptr z, gr_ctx_t ctx) { return GR_TERNARY_OP(ctx, LAGUERRE_L)(res, n, m, z, ctx); }
GR_SPECIAL_DEF int gr_hermite_h(gr_ptr res, gr_srcptr n, gr_srcptr z, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, HERMITE_H)(res, n, z, ctx); }
GR_SPECIAL_DEF int gr_legendre_p(gr_ptr res, gr_srcptr n, gr_srcptr m, gr_srcptr z, int type, gr_ctx_t ctx) { return GR_TERNARY_OP_WITH_FLAG(ctx, LEGENDRE_P)(res, n, m, z, type, ctx); }
GR_SPECIAL_DEF int gr_legendre_q(gr_ptr res, gr_srcptr n, gr_srcptr m, gr_srcptr z, int type, gr_ctx_t ctx) { return GR_TERNARY_OP_WITH_FLAG(ctx, LEGENDRE_Q)(res, n, m, z, type, ctx); }
GR_SPECIAL_DEF int gr_spherical_y_si(gr_ptr res, slong n, slong m, gr_srcptr theta, gr_srcptr phi, gr_ctx_t ctx) { return GR_SI_SI_QUATERNARY_OP(ctx, SPHERICAL_Y_SI)(res, n, m, theta, phi, ctx); }
GR_SPECIAL_DEF int gr_legendre_p_root_ui(gr_ptr root, gr_ptr weight, ulong n, ulong k, gr_ctx_t ctx) { return GR_BINARY_BINARY_OP_UI_UI(ctx, LEGENDRE_P_ROOT_UI)(root, weight, n, k, ctx); }

typedef int ((*gr_method_pfq_op_op)(gr_ptr, const gr_vec_t, const gr_vec_t, gr_srcptr, int, gr_ctx_ptr));
#define GR_PFQ_OP(ctx, NAME) (((gr_method_pfq_op_op *) ctx->methods)[GR_METHOD_ ## NAME])

GR_SPECIAL_DEF int gr_hypgeom_0f1(gr_ptr res, gr_srcptr a, gr_srcptr z, int flags, gr_ctx_t ctx) { return GR_BINARY_OP_WITH_FLAG(ctx, HYPGEOM_0F1)(res, a, z, flags, ctx); }
GR_SPECIAL_DEF int gr_hypgeom_1f1(gr_ptr res, gr_srcptr a, gr_srcptr b, gr_srcptr z, int flags, gr_ctx_t ctx) { return GR_TERNARY_OP_WITH_FLAG(ctx, HYPGEOM_1F1)(res, a, b, z, flags, ctx); }
GR_SPECIAL_DEF int gr_hypgeom_u(gr_ptr res, gr_srcptr a, gr_srcptr b, gr_srcptr z, int flags, gr_ctx_t ctx) { return GR_TERNARY_OP_WITH_FLAG(ctx, HYPGEOM_U)(res, a, b, z, flags, ctx); }
GR_SPECIAL_DEF int gr_hypgeom_2f1(gr_ptr res, gr_srcptr a, gr_srcptr b, gr_srcptr c, gr_srcptr z, int flags, gr_ctx_t ctx) { return GR_QUATERNARY_OP_WITH_FLAG(ctx, HYPGEOM_2F1)(res, a, b, c, z, flags, ctx); }
GR_SPECIAL_DEF int gr_hypgeom_pfq(gr_ptr res, const gr_vec_t a, const gr_vec_t b, gr_srcptr z, int flags, gr_ctx_t ctx) { return GR_PFQ_OP(ctx, HYPGEOM_PFQ)(res, a, b, z, flags, ctx); }

GR_SPECIAL_DEF int gr_zeta(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ZETA)(res, x, ctx); }
GR_SPECIAL_DEF int gr_zeta_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_UI(ctx, ZETA_UI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_hurwitz_zeta(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, HURWITZ_ZETA)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_polylog(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, POLYLOG)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_polygamma(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, POLYGAMMA)(res, x, y, ctx); }
GR_SPECIAL_DEF int gr_lerch_phi(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_ctx_t ctx) { return GR_TERNARY_OP(ctx, LERCH_PHI)(res, x, y, z, ctx); }
GR_SPECIAL_DEF int gr_stieltjes(gr_ptr res, const fmpz_t x, gr_srcptr y, gr_ctx_t ctx) { return GR_FMPZ_BINARY_OP(ctx, STIELTJES)(res, x, y, ctx); }

GR_SPECIAL_DEF int gr_dirichlet_eta(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, DIRICHLET_ETA)(res, x, ctx); }
GR_SPECIAL_DEF int gr_dirichlet_beta(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, DIRICHLET_BETA)(res, x, ctx); }
GR_SPECIAL_DEF int gr_riemann_xi(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, RIEMANN_XI)(res, x, ctx); }
GR_SPECIAL_DEF int gr_zeta_zero(gr_ptr res, const fmpz_t n, gr_ctx_t ctx) { return GR_UNARY_OP_FMPZ(ctx, ZETA_ZERO)(res, n, ctx); }
GR_SPECIAL_DEF int gr_zeta_zero_vec(gr_ptr res, const fmpz_t n, slong len, gr_ctx_t ctx) { return GR_BINARY_OP_FMPZ_SI(ctx, ZETA_ZERO_VEC)(res, n, len, ctx); }
GR_SPECIAL_DEF int gr_zeta_nzeros(gr_ptr res, gr_srcptr t, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ZETA_NZEROS)(res, t, ctx); }

/* todo: backlund_s, etc */

#ifdef DIRICHLET_H

WARN_UNUSED_RESULT int gr_dirichlet_chi_fmpz(gr_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, const fmpz_t n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_dirichlet_chi_vec(gr_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, gr_ctx_t ctx);  /* todo: wrap, test */
WARN_UNUSED_RESULT int gr_dirichlet_l(gr_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, gr_srcptr s, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_dirichlet_l_all(gr_vec_t res, const dirichlet_group_t G, gr_srcptr s, gr_ctx_t ctx);  /* todo: wrap, test */
WARN_UNUSED_RESULT int gr_dirichlet_hardy_theta(gr_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, gr_srcptr t, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_dirichlet_hardy_z(gr_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, gr_srcptr t, gr_ctx_t ctx);

#endif

GR_SPECIAL_DEF int gr_jacobi_theta(gr_ptr res1, gr_ptr res2, gr_ptr res3, gr_ptr res4, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx) { return GR_QUATERNARY_BINARY_OP(ctx, JACOBI_THETA)(res1, res2, res3, res4, z, tau, ctx); }
GR_SPECIAL_DEF int gr_jacobi_theta_1(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, JACOBI_THETA_1)(res, z, tau, ctx); }
GR_SPECIAL_DEF int gr_jacobi_theta_2(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, JACOBI_THETA_2)(res, z, tau, ctx); }
GR_SPECIAL_DEF int gr_jacobi_theta_3(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, JACOBI_THETA_3)(res, z, tau, ctx); }
GR_SPECIAL_DEF int gr_jacobi_theta_4(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, JACOBI_THETA_4)(res, z, tau, ctx); }

/*
GR_SPECIAL_DEF int gr_jacobi_theta_q(gr_ptr res1, gr_ptr res2, gr_ptr res3, gr_ptr res4, gr_srcptr w, gr_srcptr q, gr_ctx_t ctx) { return GR_QUATERNARY_BINARY_OP(ctx, JACOBI_THETA_Q)(res1, res2, res3, res4, w, q, ctx); }
GR_SPECIAL_DEF int gr_jacobi_theta_q_1(gr_ptr res, gr_srcptr w, gr_srcptr q, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, JACOBI_THETA_Q_1)(res, w, q, ctx); }
GR_SPECIAL_DEF int gr_jacobi_theta_q_2(gr_ptr res, gr_srcptr w, gr_srcptr q, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, JACOBI_THETA_Q_2)(res, w, q, ctx); }
GR_SPECIAL_DEF int gr_jacobi_theta_q_3(gr_ptr res, gr_srcptr w, gr_srcptr q, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, JACOBI_THETA_Q_3)(res, w, q, ctx); }
GR_SPECIAL_DEF int gr_jacobi_theta_q_4(gr_ptr res, gr_srcptr w, gr_srcptr q, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, JACOBI_THETA_Q_4)(res, w, q, ctx); }
*/

GR_SPECIAL_DEF int gr_modular_j(gr_ptr res, gr_srcptr tau, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, MODULAR_J)(res, tau, ctx); }
GR_SPECIAL_DEF int gr_modular_lambda(gr_ptr res, gr_srcptr tau, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, MODULAR_LAMBDA)(res, tau, ctx); }
GR_SPECIAL_DEF int gr_modular_delta(gr_ptr res, gr_srcptr tau, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, MODULAR_DELTA)(res, tau, ctx); }

GR_SPECIAL_DEF int gr_hilbert_class_poly(gr_ptr res, slong D, gr_srcptr x, gr_ctx_t ctx) { return GR_SI_BINARY_OP(ctx, HILBERT_CLASS_POLY)(res, D, x, ctx); }

GR_SPECIAL_DEF int gr_dedekind_eta(gr_ptr res, gr_srcptr tau, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, DEDEKIND_ETA)(res, tau, ctx); }
GR_SPECIAL_DEF int gr_dedekind_eta_q(gr_ptr res, gr_srcptr tau, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, DEDEKIND_ETA_Q)(res, tau, ctx); }

GR_SPECIAL_DEF int gr_eisenstein_e(gr_ptr res, ulong n, gr_srcptr tau, gr_ctx_t ctx) { return GR_UI_BINARY_OP(ctx, EISENSTEIN_E)(res, n, tau, ctx); }
GR_SPECIAL_DEF int gr_eisenstein_g(gr_ptr res, ulong n, gr_srcptr tau, gr_ctx_t ctx) { return GR_UI_BINARY_OP(ctx, EISENSTEIN_G)(res, n, tau, ctx); }
GR_SPECIAL_DEF int gr_eisenstein_g_vec(gr_ptr res, gr_srcptr tau, slong len, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, EISENSTEIN_G_VEC)(res, tau, len, ctx); }

GR_SPECIAL_DEF int gr_agm1(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, AGM1)(res, x, ctx); }
GR_SPECIAL_DEF int gr_agm(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, AGM)(res, x, y, ctx); }

GR_SPECIAL_DEF int gr_elliptic_k(gr_ptr res, gr_srcptr m, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ELLIPTIC_K)(res, m, ctx); }
GR_SPECIAL_DEF int gr_elliptic_e(gr_ptr res, gr_srcptr m, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ELLIPTIC_E)(res, m, ctx); }
GR_SPECIAL_DEF int gr_elliptic_pi(gr_ptr res, gr_srcptr n, gr_srcptr m, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, ELLIPTIC_PI)(res, n, m, ctx); }
GR_SPECIAL_DEF int gr_elliptic_f(gr_ptr res, gr_srcptr phi, gr_srcptr m, int pi, gr_ctx_t ctx) { return GR_BINARY_OP_WITH_FLAG(ctx, ELLIPTIC_F)(res, phi, m, pi, ctx); }
GR_SPECIAL_DEF int gr_elliptic_e_inc(gr_ptr res, gr_srcptr phi, gr_srcptr m, int pi, gr_ctx_t ctx) { return GR_BINARY_OP_WITH_FLAG(ctx, ELLIPTIC_E_INC)(res, phi, m, pi, ctx); }
GR_SPECIAL_DEF int gr_elliptic_pi_inc(gr_ptr res, gr_srcptr n, gr_srcptr phi, gr_srcptr m, int pi, gr_ctx_t ctx) { return GR_TERNARY_OP_WITH_FLAG(ctx, ELLIPTIC_PI_INC)(res, n, phi, m, pi, ctx); }

GR_SPECIAL_DEF int gr_carlson_rc(gr_ptr res, gr_srcptr x, gr_srcptr y, int flags, gr_ctx_t ctx) { return GR_BINARY_OP_WITH_FLAG(ctx, CARLSON_RC)(res, x, y, flags, ctx); }
GR_SPECIAL_DEF int gr_carlson_rf(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, int flags, gr_ctx_t ctx) { return GR_TERNARY_OP_WITH_FLAG(ctx, CARLSON_RF)(res, x, y, z, flags, ctx); }
GR_SPECIAL_DEF int gr_carlson_rd(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, int flags, gr_ctx_t ctx) { return GR_TERNARY_OP_WITH_FLAG(ctx, CARLSON_RD)(res, x, y, z, flags, ctx); }
GR_SPECIAL_DEF int gr_carlson_rg(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, int flags, gr_ctx_t ctx) { return GR_TERNARY_OP_WITH_FLAG(ctx, CARLSON_RG)(res, x, y, z, flags, ctx); }
GR_SPECIAL_DEF int gr_carlson_rj(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_srcptr z, gr_srcptr w, int flags, gr_ctx_t ctx) { return GR_QUATERNARY_OP_WITH_FLAG(ctx, CARLSON_RJ)(res, x, y, z, w, flags, ctx); }

GR_SPECIAL_DEF int gr_elliptic_invariants(gr_ptr res1, gr_ptr res2, gr_srcptr tau, gr_ctx_t ctx) { return GR_BINARY_UNARY_OP(ctx, ELLIPTIC_INVARIANTS)(res1, res2, tau, ctx); }
GR_SPECIAL_DEF int gr_elliptic_roots(gr_ptr res1, gr_ptr res2, gr_ptr res3, gr_srcptr tau, gr_ctx_t ctx) { return GR_TERNARY_UNARY_OP(ctx, ELLIPTIC_ROOTS)(res1, res2, res3, tau, ctx); }

GR_SPECIAL_DEF int gr_weierstrass_p(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, WEIERSTRASS_P)(res, z, tau, ctx); }
GR_SPECIAL_DEF int gr_weierstrass_p_prime(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, WEIERSTRASS_P_PRIME)(res, z, tau, ctx); }
GR_SPECIAL_DEF int gr_weierstrass_p_inv(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, WEIERSTRASS_P_INV)(res, z, tau, ctx); }
GR_SPECIAL_DEF int gr_weierstrass_zeta(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, WEIERSTRASS_ZETA)(res, z, tau, ctx); }
GR_SPECIAL_DEF int gr_weierstrass_sigma(gr_ptr res, gr_srcptr z, gr_srcptr tau, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, WEIERSTRASS_SIGMA)(res, z, tau, ctx); }


/* todo


    bernoulli
    eulerpoly
    zeta_ui_vec

    cyclotomic
    dedekind_sum
    moebius, ...
    prime_pi, primorial
    etc


*/

/* generic implementations */

WARN_UNUSED_RESULT int gr_generic_exp(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_expm1(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_exp2(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_exp10(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_log(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_log1p(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_log2(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_log10(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_sin(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_cos(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_sin_cos(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_tan(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_asin(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_asinh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_atan(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_atanh(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_acot(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_asec(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_acsc(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_acoth(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_asech(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_acsch(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_fac(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_fac_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_fac_ui(gr_ptr res, ulong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_fac_vec(gr_ptr res, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_rfac(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_rfac_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_rfac_ui(gr_ptr res, ulong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_rfac_vec(gr_ptr res, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_rising(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_rising_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_falling(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_falling_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_bin(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_bin_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_bin_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_bin_vec(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_bin_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_doublefac(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_doublefac_ui(gr_ptr res, ulong n, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_harmonic(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_harmonic_ui(gr_ptr res, ulong n, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_beta(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_fib2_fmpz(gr_ptr v, gr_ptr u, const fmpz_t n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_fib_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_fib_ui(gr_ptr res, ulong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_fib_vec(gr_ptr res, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_bellnum_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_bellnum_ui(gr_ptr res, ulong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_bellnum_vec(gr_ptr res, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_partitions_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_partitions_ui(gr_ptr res, ulong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_partitions_vec(gr_ptr res, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_chebyshev_t2_fmpz(gr_ptr a, gr_ptr b, const fmpz_t n, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_chebyshev_t_fmpz(gr_ptr y, const fmpz_t n, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_chebyshev_u2_fmpz(gr_ptr a, gr_ptr b, const fmpz_t n, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_chebyshev_u_fmpz(gr_ptr y, const fmpz_t n, gr_srcptr x, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_erfcx(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_hilbert_class_poly(gr_ptr res, slong D, gr_srcptr x, gr_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
