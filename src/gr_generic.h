/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_GENERIC_H
#define GR_GENERIC_H

#ifdef GR_GENERIC_INLINES_C
#define GR_GENERIC_INLINE
#else
#define GR_GENERIC_INLINE static inline
#endif

#include "flint.h"
#include "gr.h"

#ifdef __cplusplus
 extern "C" {
#endif

#define GR_GENERIC_DEBUG_RINGS 0

#if GR_GENERIC_DEBUG_RINGS
void gr_generic_init(void) { flint_throw(FLINT_ERROR, "ctx must implement init()\n"); }
void gr_generic_clear(void) { flint_throw(FLINT_ERROR, "ctx must implement clear()\n"); }
void gr_generic_swap(void) { flint_throw(FLINT_ERROR, "ctx must implement swap()\n"); }
void gr_generic_randtest(void) { flint_throw(FLINT_ERROR, "ctx must implement randtest()\n"); }
void gr_generic_write(void) { flint_throw(FLINT_ERROR, "ctx must implement write()\n"); }
void gr_generic_zero(void) { flint_throw(FLINT_ERROR, "ctx must implement zero()\n"); }
void gr_generic_one(void) { flint_throw(FLINT_ERROR, "ctx must implement one()\n"); }
void gr_generic_equal(void) { flint_throw(FLINT_ERROR, "ctx must implement equal()\n"); }
void gr_generic_set(void) { flint_throw(FLINT_ERROR, "ctx must implement set()\n"); }
void gr_generic_set_si(void) { flint_throw(FLINT_ERROR, "ctx must implement set_si()\n"); }
void gr_generic_set_ui(void) { flint_throw(FLINT_ERROR, "ctx must implement set_ui()\n"); }
void gr_generic_set_fmpz(void) { flint_throw(FLINT_ERROR, "ctx must implement set_fmpz()\n"); }
void gr_generic_neg(void) { flint_throw(FLINT_ERROR, "ctx must implement neg()\n"); }
void gr_generic_add(void) { flint_throw(FLINT_ERROR, "ctx must implement add()\n"); }
void gr_generic_sub(void) { flint_throw(FLINT_ERROR, "ctx must implement sub()\n"); }
void gr_generic_mul(void) { flint_throw(FLINT_ERROR, "ctx must implement mul()\n"); }
#else
#define gr_generic_init gr_not_implemented
#define gr_generic_clear gr_not_implemented
#define gr_generic_swap gr_not_implemented
#define gr_generic_randtest gr_not_implemented
#define gr_generic_write gr_not_implemented
#define gr_generic_zero gr_not_implemented
#define gr_generic_one gr_not_implemented
#define gr_generic_equal gr_not_implemented
#define gr_generic_set gr_not_implemented
#define gr_generic_set_si gr_not_implemented
#define gr_generic_set_ui gr_not_implemented
#define gr_generic_set_fmpz gr_not_implemented
#define gr_generic_neg gr_not_implemented
#define gr_generic_add gr_not_implemented
#define gr_generic_sub gr_not_implemented
#define gr_generic_mul gr_not_implemented
#endif



/* some useful generic functions, currently not overloadable */
#ifdef FMPZ_POLY_H
WARN_UNUSED_RESULT int _gr_fmpz_poly_evaluate_horner(gr_ptr res, const fmpz * f, slong len, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_fmpz_poly_evaluate_horner(gr_ptr res, const fmpz_poly_t f, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_fmpz_poly_evaluate_rectangular(gr_ptr res, const fmpz * f, slong len, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_fmpz_poly_evaluate_rectangular(gr_ptr res, const fmpz_poly_t f, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_fmpz_poly_evaluate(gr_ptr res, const fmpz * f, slong len, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_fmpz_poly_evaluate(gr_ptr res, const fmpz_poly_t f, gr_srcptr x, gr_ctx_t ctx);
#endif

#ifdef FMPZ_MPOLY_H
WARN_UNUSED_RESULT int gr_fmpz_mpoly_evaluate(gr_ptr res, const fmpz_mpoly_t f, gr_srcptr x, const fmpz_mpoly_ctx_t mctx, gr_ctx_t ctx);
#endif

truth_t gr_generic_ctx_predicate(gr_ctx_t ctx);
truth_t gr_generic_ctx_predicate_true(gr_ctx_t ctx);
truth_t gr_generic_ctx_predicate_false(gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_ctx_clear(gr_ctx_t ctx);

void gr_generic_set_shallow(gr_ptr res, gr_srcptr x, const gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_write_n(gr_stream_t out, gr_srcptr x, slong n, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_randtest_not_zero(gr_ptr x, flint_rand_t state, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_randtest_small(gr_ptr x, flint_rand_t state, gr_ctx_t ctx);

WARN_UNUSED_RESULT truth_t gr_generic_is_zero(gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT truth_t gr_generic_is_one(gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT truth_t gr_generic_is_neg_one(gr_srcptr x, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_neg_one(gr_ptr res, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_set_other(gr_ptr res, gr_srcptr x, gr_ctx_t xctx, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_set_fmpq(gr_ptr res, const fmpq_t y, gr_ctx_t ctx);

#ifdef FEXPR_H
WARN_UNUSED_RESULT int gr_generic_set_fexpr(gr_ptr res, fexpr_vec_t inputs, gr_vec_t outputs, const fexpr_t expr, gr_ctx_t ctx);
#endif

WARN_UNUSED_RESULT int gr_generic_add_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_add_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_add_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_add_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_add_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_other_add(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_sub_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_sub_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_sub_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_sub_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_sub_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_other_sub(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_mul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_mul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_mul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_mul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_mul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_other_mul(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_addmul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_addmul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_addmul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_addmul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_addmul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_addmul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_submul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_submul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_submul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_submul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_submul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_submul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_mul_two(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_sqr(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_mul_2exp_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_mul_2exp_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_set_fmpz_2exp_fmpz(gr_ptr res, const fmpz_t x, const fmpz_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_get_fmpz_2exp_fmpz(fmpz_t res1, fmpz_t res2, gr_ptr x, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_inv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT truth_t gr_generic_is_invertible(gr_srcptr x, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_div_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_div_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_div_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_div_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_div_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_other_div(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_divexact(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_pow_fmpz_sliding(gr_ptr f, gr_srcptr g, const fmpz_t pow, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_pow_ui_sliding(gr_ptr f, gr_srcptr g, ulong pow, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_pow_fmpz_binexp(gr_ptr res, gr_srcptr x, const fmpz_t exp, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_pow_ui_binexp(gr_ptr res, gr_srcptr x, ulong e, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_pow_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t e, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_pow_si(gr_ptr res, gr_srcptr x, slong e, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_pow_ui(gr_ptr res, gr_srcptr x, ulong e, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_pow_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_pow_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_other_pow(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_numerator(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_denominator(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);

WARN_UNUSED_RESULT truth_t gr_generic_is_square(gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_sqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_rsqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_cmp(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_cmpabs(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_cmp_other(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_cmpabs_other(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_bernoulli_ui(gr_ptr res, ulong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_bernoulli_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_bernoulli_vec(gr_ptr res, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_eulernum_ui(gr_ptr res, ulong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_eulernum_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_eulernum_vec(gr_ptr res, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_stirling_s1u_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_stirling_s1_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_stirling_s2_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_stirling_s1u_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_stirling_s1_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_stirling_s2_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx);

void gr_generic_vec_init(gr_ptr vec, slong len, gr_ctx_t ctx);
void gr_generic_vec_clear(gr_ptr vec, slong len, gr_ctx_t ctx);
void gr_generic_vec_swap(gr_ptr vec1, gr_ptr vec2, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_vec_zero(gr_ptr vec, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_set(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_neg(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_normalise(slong * res, gr_srcptr vec, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT slong gr_generic_vec_normalise_weak(gr_srcptr vec, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_mul_scalar_2exp_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_scalar_addmul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_scalar_submul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_scalar_addmul_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_scalar_submul_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT truth_t gr_generic_vec_equal(gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_is_zero(gr_srcptr vec, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_dot(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_dot_rev(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_dot_ui(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const ulong * vec2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_dot_si(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const slong * vec2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_dot_fmpz(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const fmpz * vec2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_set_powers(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_reciprocals(gr_ptr res, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_pow_fmpz_sliding(gr_ptr f, gr_srcptr g, const fmpz_t pow, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_pow_ui_sliding(gr_ptr f, gr_srcptr g, ulong pow, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_pow_fmpz_binexp(gr_ptr res, gr_srcptr x, const fmpz_t exp, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_pow_ui_binexp(gr_ptr res, gr_srcptr x, ulong e, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_vec_add(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_sub(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_mul(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_div(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_divexact(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_pow(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_vec_add_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_sub_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_mul_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_div_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_divexact_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_pow_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_vec_add_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_sub_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_mul_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_div_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_divexact_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_pow_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_vec_add_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_sub_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_mul_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_div_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_divexact_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_pow_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_vec_add_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_sub_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_mul_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_div_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_divexact_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_pow_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_vec_add_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_sub_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_mul_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_div_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_divexact_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_pow_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_scalar_add_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_scalar_sub_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_scalar_mul_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_scalar_div_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_scalar_divexact_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_scalar_pow_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_vec_add_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_sub_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_mul_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_div_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_divexact_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_pow_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_other_add_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_other_sub_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_other_mul_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_other_div_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_other_divexact_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_other_pow_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_vec_add_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_sub_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_mul_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_div_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_divexact_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_vec_pow_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_generic_scalar_other_add_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_scalar_other_sub_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_scalar_other_mul_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_scalar_other_div_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_scalar_other_divexact_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_generic_scalar_other_pow_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx);


#ifdef __cplusplus
}
#endif

#endif
