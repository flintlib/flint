/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef GR_POLY_H
#define GR_POLY_H

#ifdef GR_POLY_INLINES_C
#define GR_POLY_INLINE FLINT_DLL
#else
#define GR_POLY_INLINE static __inline__
#endif

#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "gr.h"

#ifdef __cplusplus
 extern "C" {
#endif

/* fixme: compatible with flint polys but not with arb, ... */
typedef struct
{
    gr_ptr coeffs;
    slong alloc;
    slong length;
}
gr_poly_struct;

typedef gr_poly_struct gr_poly_t[1];

void gr_poly_init(gr_poly_t poly, gr_ctx_t ctx);
void gr_poly_init2(gr_poly_t poly, slong len, gr_ctx_t ctx);
void gr_poly_clear(gr_poly_t poly, gr_ctx_t ctx);

GR_POLY_INLINE gr_ptr
gr_poly_entry_ptr(gr_poly_t poly, slong i, gr_ctx_t ctx)
{
    return GR_ENTRY(poly->coeffs, i, ctx->sizeof_elem);
}

GR_POLY_INLINE slong gr_poly_length(const gr_poly_t poly, gr_ctx_t ctx)
{
    return poly->length;
}

GR_POLY_INLINE void
gr_poly_swap(gr_poly_t poly1, gr_poly_t poly2, gr_ctx_t ctx)
{
    gr_poly_struct t = *poly1;
    *poly1 = *poly2;
    *poly2 = t;
}

void gr_poly_fit_length(gr_poly_t poly, slong len, gr_ctx_t ctx);
void _gr_poly_set_length(gr_poly_t poly, slong len, gr_ctx_t ctx);
void _gr_poly_normalise(gr_poly_t poly, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_poly_set(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_reverse(gr_ptr res, gr_srcptr poly, slong len, slong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_reverse(gr_poly_t res, const gr_poly_t poly, slong n, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_poly_truncate(gr_poly_t poly, slong newlen, gr_ctx_t ctx);

GR_POLY_INLINE WARN_UNUSED_RESULT int
gr_poly_zero(gr_poly_t poly, gr_ctx_t ctx)
{
    _gr_poly_set_length(poly, 0, ctx);
    return GR_SUCCESS;
}

WARN_UNUSED_RESULT int gr_poly_one(gr_poly_t poly, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_neg_one(gr_poly_t poly, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_gen(gr_poly_t poly, gr_ctx_t ctx);

int gr_poly_write(gr_stream_t out, const gr_poly_t poly, gr_ctx_t ctx);
int gr_poly_print(const gr_poly_t poly, gr_ctx_t ctx);

int gr_poly_randtest(gr_poly_t poly, flint_rand_t state, slong len, gr_ctx_t ctx);

truth_t _gr_poly_equal(gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx);
truth_t gr_poly_equal(const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx);

truth_t gr_poly_is_zero(const gr_poly_t poly, gr_ctx_t ctx);
truth_t gr_poly_is_one(const gr_poly_t poly, gr_ctx_t ctx);
truth_t gr_poly_is_gen(const gr_poly_t poly, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_poly_set_scalar(gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_set_si(gr_poly_t poly, slong x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_set_ui(gr_poly_t poly, slong x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_set_fmpz(gr_poly_t poly, const fmpz_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_set_fmpq(gr_poly_t poly, const fmpq_t x, gr_ctx_t ctx);

/* todo: test */
WARN_UNUSED_RESULT int gr_poly_set_fmpz_poly(gr_poly_t res, const fmpz_poly_t src, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_set_fmpq_poly(gr_poly_t res, const fmpq_poly_t src, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_set_gr_poly_other(gr_poly_t res, const gr_poly_t x, gr_ctx_t x_ctx, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_poly_set_coeff_scalar(gr_poly_t poly, slong n, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_set_coeff_si(gr_poly_t poly, slong n, slong x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_set_coeff_ui(gr_poly_t poly, slong n, ulong x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_set_coeff_fmpz(gr_poly_t poly, slong n, const fmpz_t x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_set_coeff_fmpq(gr_poly_t poly, slong n, const fmpq_t x, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_poly_get_coeff_scalar(gr_ptr res, const gr_poly_t poly, slong n, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_poly_neg(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_add(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_add(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_sub(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_sub(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_mul(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_mul(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_mullow_generic(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx);
GR_POLY_INLINE WARN_UNUSED_RESULT int _gr_poly_mullow(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong len, gr_ctx_t ctx) { return GR_POLY_BINARY_TRUNC_OP(ctx, POLY_MULLOW)(res, poly1, len1, poly2, len2, len, ctx); }
WARN_UNUSED_RESULT int gr_poly_mullow(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_mul_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx);

/* powering */

WARN_UNUSED_RESULT int _gr_poly_pow_series_ui_binexp(gr_ptr res, gr_srcptr f, slong flen, ulong exp, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_pow_series_ui_binexp(gr_poly_t res, const gr_poly_t poly, ulong exp, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_pow_series_ui(gr_ptr res, gr_srcptr f, slong flen, ulong exp, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_pow_series_ui(gr_poly_t res, const gr_poly_t poly, ulong exp, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_pow_ui_binexp(gr_ptr res, gr_srcptr f, slong flen, ulong exp, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_pow_ui_binexp(gr_poly_t res, const gr_poly_t poly, ulong exp, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_pow_ui(gr_ptr res, gr_srcptr f, slong flen, ulong exp, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_pow_ui(gr_poly_t res, const gr_poly_t poly, ulong exp, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_poly_pow_fmpz(gr_poly_t res, const gr_poly_t poly, const fmpz_t exp, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_pow_series_fmpq_recurrence(gr_ptr h, gr_srcptr f, slong flen, const fmpq_t exp, slong len, int precomp, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_pow_series_fmpq_recurrence(gr_poly_t res, const gr_poly_t poly, const fmpq_t exp, slong len, gr_ctx_t ctx);

/* shifting */

WARN_UNUSED_RESULT int _gr_poly_shift_left(gr_ptr res, gr_srcptr poly, slong len, slong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_shift_left(gr_poly_t res, const gr_poly_t poly, slong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_poly_shift_right(gr_ptr res, gr_srcptr poly, slong len, slong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_shift_right(gr_poly_t res, const gr_poly_t poly, slong n, gr_ctx_t ctx);

/* division */

WARN_UNUSED_RESULT int gr_poly_inv(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_divrem_basecase_preinv(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_poly_divrem_basecase_noinv(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_poly_divrem_basecase(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_divrem_basecase(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_divrem_divconquer_preinv(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, slong cutoff, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_poly_divrem_divconquer(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong cutoff, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_divrem_divconquer(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, slong cutoff, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_div_newton(gr_ptr Q, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_poly_divrem_newton(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_divrem_newton(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx);

GR_POLY_INLINE WARN_UNUSED_RESULT int _gr_poly_divrem(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx) { return GR_POLY_BINARY_BINARY_OP(ctx, POLY_DIVREM)(Q, R, A, lenA, B, lenB, ctx); }

WARN_UNUSED_RESULT int _gr_poly_divrem_generic(gr_ptr Q, gr_ptr R, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_divrem(gr_poly_t Q, gr_poly_t R, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx);

/* todo: div with fast divisibility checking; rem; divexact */

WARN_UNUSED_RESULT int _gr_poly_inv_series_newton(gr_ptr Qinv, gr_srcptr Q, slong Qlen, slong len, slong cutoff, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_inv_series_newton(gr_poly_t Qinv, const gr_poly_t Q, slong len, slong cutoff, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_poly_inv_series_basecase(gr_ptr Qinv, gr_srcptr Q, slong Qlen, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_inv_series_basecase(gr_poly_t Qinv, const gr_poly_t Q, slong len, gr_ctx_t ctx);
GR_POLY_INLINE WARN_UNUSED_RESULT int _gr_poly_inv_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx) { return GR_POLY_UNARY_TRUNC_OP(ctx, POLY_INV_SERIES)(res, f, flen, len, ctx); }
WARN_UNUSED_RESULT int _gr_poly_inv_series_generic(gr_ptr Qinv, gr_srcptr Q, slong Qlen, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_inv_series(gr_poly_t Qinv, const gr_poly_t Q, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_div_series_newton(gr_ptr res, gr_srcptr B, slong Blen, gr_srcptr A, slong Alen, slong len, slong cutoff, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_div_series_newton(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, slong len, slong cutoff, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_poly_div_series_basecase(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_div_series_basecase(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_poly_div_series_invmul(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_div_series_invmul(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx);
GR_POLY_INLINE WARN_UNUSED_RESULT int _gr_poly_div_series(gr_ptr res, gr_srcptr f, slong flen, gr_srcptr g, slong glen, slong len, gr_ctx_t ctx) { return GR_POLY_BINARY_TRUNC_OP(ctx, POLY_DIV_SERIES)(res, f, flen, g, glen, len, ctx); }
WARN_UNUSED_RESULT int _gr_poly_div_series_generic(gr_ptr Q, gr_srcptr A, slong Alen, gr_srcptr B, slong Blen, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_div_series(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_sqrt_series_newton(gr_ptr res, gr_srcptr f, slong flen, slong len, slong cutoff, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_sqrt_series_newton(gr_poly_t res, const gr_poly_t f, slong len, slong cutoff, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_poly_sqrt_series_basecase(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_sqrt_series_basecase(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_poly_sqrt_series_miller(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_sqrt_series_miller(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx);
GR_POLY_INLINE WARN_UNUSED_RESULT int _gr_poly_sqrt_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx) { return GR_POLY_UNARY_TRUNC_OP(ctx, POLY_SQRT_SERIES)(res, f, flen, len, ctx); }
WARN_UNUSED_RESULT int _gr_poly_sqrt_series_generic(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_sqrt_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_rsqrt_series_newton(gr_ptr res, gr_srcptr f, slong flen, slong len, slong cutoff, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_rsqrt_series_newton(gr_poly_t res, const gr_poly_t f, slong len, slong cutoff, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_poly_rsqrt_series_basecase(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_rsqrt_series_basecase(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_poly_rsqrt_series_miller(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_rsqrt_series_miller(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx);
GR_POLY_INLINE WARN_UNUSED_RESULT int _gr_poly_rsqrt_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx) { return GR_POLY_UNARY_TRUNC_OP(ctx, POLY_RSQRT_SERIES)(res, f, flen, len, ctx); }
WARN_UNUSED_RESULT int _gr_poly_rsqrt_series_generic(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_rsqrt_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_evaluate_rectangular(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_evaluate_rectangular(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_evaluate_horner(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_evaluate_horner(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_evaluate(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_evaluate(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_evaluate_other_horner(gr_ptr res, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_evaluate_other_horner(gr_ptr res, const gr_poly_t f, gr_srcptr a, gr_ctx_t a_ctx, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_poly_evaluate_other_rectangular(gr_ptr res, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_evaluate_other_rectangular(gr_ptr res, const gr_poly_t f, gr_srcptr a, gr_ctx_t a_ctx, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_poly_evaluate_other(gr_ptr res, gr_srcptr f, slong len, const gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_evaluate_other(gr_ptr res, const gr_poly_t f, gr_srcptr a, gr_ctx_t a_ctx, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_taylor_shift_horner(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_taylor_shift_horner(gr_poly_t res, const gr_poly_t f, gr_srcptr c, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_taylor_shift_divconquer(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_taylor_shift_divconquer(gr_poly_t res, const gr_poly_t f, gr_srcptr c, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_taylor_shift(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_taylor_shift(gr_poly_t res, const gr_poly_t f, gr_srcptr c, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_compose_horner(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_compose_horner(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_compose_divconquer(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_compose_divconquer(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_compose(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_compose(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_derivative(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_derivative(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_integral(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_integral(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_make_monic(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_make_monic(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx);

WARN_UNUSED_RESULT truth_t _gr_poly_is_monic(gr_srcptr poly, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT truth_t gr_poly_is_monic(const gr_poly_t res, gr_ctx_t ctx);

/* GCD, resultant */

WARN_UNUSED_RESULT int _gr_poly_hgcd(slong * sgn, gr_ptr * M, slong * lenM, gr_ptr A, slong * lenA, gr_ptr B, slong * lenB, gr_srcptr a, slong lena, gr_srcptr b, slong lenb, slong cutoff, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_gcd_hgcd(gr_ptr G, slong * _lenG, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong inner_cutoff, slong cutoff, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_gcd_hgcd(gr_poly_t G, const gr_poly_t A, const gr_poly_t B, slong inner_cutoff, slong cutoff, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_gcd_euclidean(gr_ptr G, slong * lenG, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_gcd_euclidean(gr_poly_t G, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_poly_gcd(gr_ptr G, slong * lenG, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_gcd(gr_poly_t G, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_xgcd_euclidean(slong * lenG, gr_ptr G, gr_ptr S, gr_ptr T, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, gr_srcptr invB, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_xgcd_euclidean(gr_poly_t G, gr_poly_t S, gr_poly_t T, const gr_poly_t A, const gr_poly_t B, gr_ctx_t ctx);

/* todo: undocumented; deduplicate */
int _gr_poly_hgcd_res(gr_ptr r, slong * sgn, gr_ptr * M, slong * lenM,
                               gr_ptr A, slong * lenA,
                               gr_ptr B, slong * lenB,
                               gr_srcptr a, slong lena,
                               gr_srcptr b, slong lenb,
                               slong cutoff,
                               gr_ctx_t ctx);

int _gr_poly_resultant_euclidean(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx);
int gr_poly_resultant_euclidean(gr_ptr r, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx);

int _gr_poly_resultant_hgcd(gr_ptr res, gr_srcptr A, slong lenA, gr_srcptr B, slong lenB, slong inner_cutoff, slong cutoff, gr_ctx_t ctx);
int gr_poly_resultant_hgcd(gr_ptr r, const gr_poly_t f, const gr_poly_t g, slong inner_cutoff, slong cutoff, gr_ctx_t ctx);

int _gr_poly_resultant_sylvester(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx);
int gr_poly_resultant_sylvester(gr_ptr r, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx);

int _gr_poly_resultant_small(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx);
int gr_poly_resultant_small(gr_ptr r, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx);

int _gr_poly_resultant(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx);
int gr_poly_resultant(gr_ptr r, const gr_poly_t f, const gr_poly_t g, gr_ctx_t ctx);

/* Multipoint evaluation/interpolation */

gr_ptr * _gr_poly_tree_alloc(slong len, gr_ctx_t ctx);
void _gr_poly_tree_free(gr_ptr * tree, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_poly_tree_build(gr_ptr * tree, gr_srcptr roots, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_evaluate_vec_fast_precomp(gr_ptr vs, gr_srcptr poly, slong plen, gr_ptr * tree, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_evaluate_vec_fast(gr_ptr ys, gr_srcptr poly, slong plen, gr_srcptr xs, slong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_evaluate_vec_fast(gr_vec_t ys, const gr_poly_t poly, const gr_vec_t xs, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_evaluate_vec_iter(gr_ptr ys, gr_srcptr poly, slong plen, gr_srcptr xs, slong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_evaluate_vec_iter(gr_vec_t ys, const gr_poly_t poly, const gr_vec_t xs, gr_ctx_t ctx);

/* Squarefree factorization */

int gr_poly_factor_squarefree(gr_ptr c, gr_vec_t fac, gr_vec_t exp, const gr_poly_t F, gr_ctx_t ctx);
int gr_poly_squarefree_part(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx);

/* Roots */

typedef int ((*gr_poly_roots_op)(gr_vec_t, gr_vec_t, const gr_poly_t, int, gr_ctx_ptr));
#define GR_POLY_ROOTS_OP(ctx, NAME) (((gr_poly_roots_op *) ctx->methods)[GR_METHOD_ ## NAME])
GR_POLY_INLINE WARN_UNUSED_RESULT int gr_poly_roots(gr_vec_t roots, gr_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx) { return GR_POLY_ROOTS_OP(ctx, POLY_ROOTS)(roots, mult, poly, flags, ctx); }

typedef int ((*gr_poly_roots_op_other)(gr_vec_t, gr_vec_t, const gr_poly_t, gr_ctx_ptr, int, gr_ctx_ptr));
#define GR_POLY_ROOTS_OP_OTHER(ctx, NAME) (((gr_poly_roots_op_other *) ctx->methods)[GR_METHOD_ ## NAME])
GR_POLY_INLINE WARN_UNUSED_RESULT int gr_poly_roots_other(gr_vec_t roots, gr_vec_t mult, const gr_poly_t poly, gr_ctx_t poly_ctx, int flags, gr_ctx_t ctx) { return GR_POLY_ROOTS_OP_OTHER(ctx, POLY_ROOTS_OTHER)(roots, mult, poly, poly_ctx, flags, ctx); }

WARN_UNUSED_RESULT int _gr_poly_atan_series(gr_ptr res, gr_srcptr A, slong Alen, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_atan_series(gr_poly_t res, const gr_poly_t A, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_poly_atanh_series(gr_ptr res, gr_srcptr A, slong Alen, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_atanh_series(gr_poly_t res, const gr_poly_t A, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_log_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_log_series(gr_poly_t res, const gr_poly_t f, slong len, gr_ctx_t ctx);

WARN_UNUSED_RESULT int _gr_poly_exp_series_basecase(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_exp_series_basecase(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_poly_exp_series_basecase_mul(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_exp_series_basecase_mul(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_poly_exp_series_newton(gr_ptr f, gr_ptr g, gr_srcptr h, slong hlen, slong n, slong cutoff, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_exp_series_newton(gr_poly_t f, const gr_poly_t h, slong n, slong cutoff, gr_ctx_t ctx);
GR_POLY_INLINE WARN_UNUSED_RESULT int _gr_poly_exp_series(gr_ptr res, gr_srcptr f, slong flen, slong len, gr_ctx_t ctx) { return GR_POLY_UNARY_TRUNC_OP(ctx, POLY_EXP_SERIES)(res, f, flen, len, ctx); }
WARN_UNUSED_RESULT int _gr_poly_exp_series_generic(gr_ptr f, gr_srcptr h, slong hlen, slong n, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_poly_exp_series(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
