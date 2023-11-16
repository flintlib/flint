/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARB_POLY_H
#define ARB_POLY_H

#ifdef ARB_POLY_INLINES_C
#define ARB_POLY_INLINE
#else
#define ARB_POLY_INLINE static inline
#endif

#include "fmpq_types.h"
#include "arb.h"
#include "acb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Memory management */

void arb_poly_init(arb_poly_t poly);

void arb_poly_init2(arb_poly_t poly, slong len);

void arb_poly_clear(arb_poly_t poly);

void arb_poly_fit_length(arb_poly_t poly, slong len);

void _arb_poly_set_length(arb_poly_t poly, slong len);

void _arb_poly_normalise(arb_poly_t poly);

ARB_POLY_INLINE void
arb_poly_swap(arb_poly_t poly1, arb_poly_t poly2)
{
    arb_poly_struct t = *poly1;
    *poly1 = *poly2;
    *poly2 = t;
}

void arb_poly_set(arb_poly_t poly, const arb_poly_t src);

void arb_poly_set_round(arb_poly_t poly, const arb_poly_t src, slong prec);

void arb_poly_set_trunc(arb_poly_t res, const arb_poly_t poly, slong n);

void arb_poly_set_trunc_round(arb_poly_t res, const arb_poly_t poly, slong n, slong prec);

/* Basic manipulation */

ARB_POLY_INLINE slong arb_poly_length(const arb_poly_t poly)
{
    return poly->length;
}

ARB_POLY_INLINE slong arb_poly_degree(const arb_poly_t poly)
{
    return poly->length - 1;
}

slong arb_poly_valuation(const arb_poly_t poly);

ARB_POLY_INLINE int
arb_poly_is_zero(const arb_poly_t z)
{
    return arb_poly_length(z) == 0;
}

ARB_POLY_INLINE int
arb_poly_is_one(const arb_poly_t z)
{
    return (arb_poly_length(z) == 1) && arb_is_one(z->coeffs);
}

ARB_POLY_INLINE int
arb_poly_is_x(const arb_poly_t z)
{
    return (arb_poly_length(z) == 2) && arb_is_zero(z->coeffs)
        && arb_is_one(z->coeffs + 1);
}

ARB_POLY_INLINE void arb_poly_zero(arb_poly_t poly)
{
    poly->length = 0;
}

ARB_POLY_INLINE void
arb_poly_one(arb_poly_t poly)
{
    arb_poly_fit_length(poly, 1);
    arb_one(poly->coeffs);
    _arb_poly_set_length(poly, 1);
}

void arb_poly_set_coeff_si(arb_poly_t poly, slong n, slong x);

void arb_poly_set_coeff_arb(arb_poly_t poly, slong n, const arb_t x);

void arb_poly_get_coeff_arb(arb_t x, const arb_poly_t poly, slong n);

#define arb_poly_get_coeff_ptr(poly, n) \
    ((n) < (poly)->length ? (poly)->coeffs + (n) : NULL)

void _arb_poly_reverse(arb_ptr res, arb_srcptr poly, slong len, slong n);

void _arb_poly_shift_right(arb_ptr res, arb_srcptr poly, slong len, slong n);

void arb_poly_shift_right(arb_poly_t res, const arb_poly_t poly, slong n);

void _arb_poly_shift_left(arb_ptr res, arb_srcptr poly, slong len, slong n);

void arb_poly_shift_left(arb_poly_t res, const arb_poly_t poly, slong n);

ARB_POLY_INLINE void
arb_poly_truncate(arb_poly_t poly, slong newlen)
{
    if (poly->length > newlen)
    {
        slong i;
        for (i = newlen; i < poly->length; i++)
            arb_zero(poly->coeffs + i);
        poly->length = newlen;
        _arb_poly_normalise(poly);
    }
}

/* Conversions */

void arb_poly_set_fmpz_poly(arb_poly_t poly, const fmpz_poly_t src, slong prec);

void arb_poly_set_fmpq_poly(arb_poly_t poly, const fmpq_poly_t src, slong prec);

ARB_POLY_INLINE void
arb_poly_set_arb(arb_poly_t poly, const arb_t c)
{
    arb_poly_fit_length(poly, 1);
    arb_set(poly->coeffs, c);
    _arb_poly_set_length(poly, !arb_is_zero(poly->coeffs));
}

void arb_poly_set_si(arb_poly_t poly, slong c);

int arb_poly_get_unique_fmpz_poly(fmpz_poly_t res, const arb_poly_t src);

/* Comparisons */

int arb_poly_contains(const arb_poly_t poly1, const arb_poly_t poly2);

int arb_poly_contains_fmpz_poly(const arb_poly_t poly1, const fmpz_poly_t poly2);

int arb_poly_contains_fmpq_poly(const arb_poly_t poly1, const fmpq_poly_t poly2);

int arb_poly_equal(const arb_poly_t A, const arb_poly_t B);

int _arb_poly_overlaps(arb_srcptr poly1, slong len1, arb_srcptr poly2, slong len2);

int arb_poly_overlaps(const arb_poly_t poly1, const arb_poly_t poly2);

/* Bounds */

void _arb_poly_majorant(arb_ptr res, arb_srcptr vec, slong len, slong prec);

void arb_poly_majorant(arb_poly_t res, const arb_poly_t poly, slong prec);

/* IO */

#ifdef FLINT_HAVE_FILE
void arb_poly_fprintd(FILE * file, const arb_poly_t poly, slong digits);
#endif

void arb_poly_printd(const arb_poly_t poly, slong digits);

/* Random generation */

void arb_poly_randtest(arb_poly_t poly, flint_rand_t state, slong len, slong prec, slong mag_bits);

/* Arithmetic */

void
_arb_poly_add(arb_ptr res, arb_srcptr poly1, slong len1,
    arb_srcptr poly2, slong len2, slong prec);

void arb_poly_add(arb_poly_t res, const arb_poly_t poly1,
              const arb_poly_t poly2, slong prec);

void arb_poly_add_si(arb_poly_t res, const arb_poly_t poly, slong c, slong prec);

void _arb_poly_sub(arb_ptr res, arb_srcptr poly1, slong len1,
    arb_srcptr poly2, slong len2, slong prec);

void arb_poly_sub(arb_poly_t res, const arb_poly_t poly1,
              const arb_poly_t poly2, slong prec);

void arb_poly_add_series(arb_poly_t res, const arb_poly_t poly1,
              const arb_poly_t poly2, slong len, slong prec);

void arb_poly_sub_series(arb_poly_t res, const arb_poly_t poly1,
              const arb_poly_t poly2, slong len, slong prec);

ARB_POLY_INLINE void
arb_poly_neg(arb_poly_t res, const arb_poly_t poly)
{
    arb_poly_fit_length(res, poly->length);
    _arb_vec_neg(res->coeffs, poly->coeffs, poly->length);
    _arb_poly_set_length(res, poly->length);
}

ARB_POLY_INLINE void
arb_poly_scalar_mul_2exp_si(arb_poly_t res, const arb_poly_t poly, slong c)
{
    arb_poly_fit_length(res, poly->length);
    _arb_vec_scalar_mul_2exp_si(res->coeffs, poly->coeffs, poly->length, c);
    _arb_poly_set_length(res, poly->length);
}

ARB_POLY_INLINE void
arb_poly_scalar_mul(arb_poly_t res, const arb_poly_t poly, const arb_t c, slong prec)
{
    arb_poly_fit_length(res, poly->length);
    _arb_vec_scalar_mul(res->coeffs, poly->coeffs, poly->length, c, prec);
    _arb_poly_set_length(res, poly->length);
    _arb_poly_normalise(res);
}

ARB_POLY_INLINE void
arb_poly_scalar_div(arb_poly_t res, const arb_poly_t poly, const arb_t c, slong prec)
{
    arb_poly_fit_length(res, poly->length);
    _arb_vec_scalar_div(res->coeffs, poly->coeffs, poly->length, c, prec);
    _arb_poly_set_length(res, poly->length);
    _arb_poly_normalise(res);
}

void _arb_poly_mullow_ztrunc(arb_ptr res,
    arb_srcptr poly1, slong len1,
    arb_srcptr poly2, slong len2, slong n, slong prec);

void arb_poly_mullow_ztrunc(arb_poly_t res, const arb_poly_t poly1,
                                            const arb_poly_t poly2,
                                                slong n, slong prec);

void _arb_poly_mullow_classical(arb_ptr res,
    arb_srcptr poly1, slong len1,
    arb_srcptr poly2, slong len2, slong n, slong prec);

void arb_poly_mullow_classical(arb_poly_t res, const arb_poly_t poly1,
                                            const arb_poly_t poly2,
                                                slong n, slong prec);

void _arb_poly_mullow_block(arb_ptr C,
    arb_srcptr A, slong lenA,
    arb_srcptr B, slong lenB, slong n, slong prec);

void arb_poly_mullow_block(arb_poly_t res, const arb_poly_t poly1,
              const arb_poly_t poly2, slong len, slong prec);

void _arb_poly_mullow(arb_ptr C,
    arb_srcptr A, slong lenA,
    arb_srcptr B, slong lenB, slong n, slong prec);

void arb_poly_mullow(arb_poly_t res, const arb_poly_t poly1,
              const arb_poly_t poly2, slong len, slong prec);

void _arb_poly_mul(arb_ptr C,
    arb_srcptr A, slong lenA,
    arb_srcptr B, slong lenB, slong prec);

void arb_poly_mul(arb_poly_t res, const arb_poly_t poly1,
              const arb_poly_t poly2, slong prec);

ARB_POLY_INLINE void
_arb_poly_mul_monic(arb_ptr res, arb_srcptr poly1, slong len1,
    arb_srcptr poly2, slong len2, slong prec)
{
    if (len1 + len2 - 2 > 0)
        _arb_poly_mullow(res, poly1, len1, poly2, len2, len1 + len2 - 2, prec);
    arb_one(res + len1 + len2 - 2);
}

void _arb_poly_inv_series(arb_ptr Qinv,
    arb_srcptr Q, slong Qlen, slong len, slong prec);

void arb_poly_inv_series(arb_poly_t Qinv, const arb_poly_t Q, slong n, slong prec);

void  _arb_poly_div_series(arb_ptr Q, arb_srcptr A, slong Alen,
    arb_srcptr B, slong Blen, slong n, slong prec);

void arb_poly_div_series(arb_poly_t Q, const arb_poly_t A, const arb_poly_t B, slong n, slong prec);

void
_arb_poly_div(arb_ptr Q,
    arb_srcptr A, slong lenA,
    arb_srcptr B, slong lenB, slong prec);

void _arb_poly_divrem(arb_ptr Q, arb_ptr R,
    arb_srcptr A, slong lenA,
    arb_srcptr B, slong lenB, slong prec);

void _arb_poly_rem(arb_ptr R,
    arb_srcptr A, slong lenA,
    arb_srcptr B, slong lenB, slong prec);

int arb_poly_divrem(arb_poly_t Q, arb_poly_t R,
                             const arb_poly_t A, const arb_poly_t B, slong prec);

void _arb_poly_div_root(arb_ptr Q, arb_t R, arb_srcptr A,
    slong len, const arb_t c, slong prec);

/* Product trees */

void _arb_poly_product_roots(arb_ptr poly, arb_srcptr xs, slong n, slong prec);

void arb_poly_product_roots(arb_poly_t poly, arb_srcptr xs, slong n, slong prec);

void _arb_poly_product_roots_complex(arb_ptr poly, arb_srcptr r, slong rn,
    acb_srcptr c, slong cn, slong prec);

void arb_poly_product_roots_complex(arb_poly_t poly,
        arb_srcptr r, slong rn, acb_srcptr c, slong cn, slong prec);

arb_ptr * _arb_poly_tree_alloc(slong len);

void _arb_poly_tree_free(arb_ptr * tree, slong len);

void _arb_poly_tree_build(arb_ptr * tree, arb_srcptr roots, slong len, slong prec);

/* Composition */

void _arb_poly_taylor_shift(arb_ptr poly, const arb_t c, slong n, slong prec);
void arb_poly_taylor_shift(arb_poly_t g, const arb_poly_t f, const arb_t c, slong prec);
void _arb_poly_compose(arb_ptr res, arb_srcptr poly1, slong len1, arb_srcptr poly2, slong len2, slong prec);
void arb_poly_compose(arb_poly_t res, const arb_poly_t poly1, const arb_poly_t poly2, slong prec);

void _arb_poly_compose_series(arb_ptr res, arb_srcptr poly1, slong len1,
                            arb_srcptr poly2, slong len2, slong n, slong prec);

void arb_poly_compose_series(arb_poly_t res,
                    const arb_poly_t poly1,
                    const arb_poly_t poly2, slong n, slong prec);

/* Reversion */

void _arb_poly_revert_series(arb_ptr Qinv, arb_srcptr Q, slong Qlen, slong n, slong prec);
void arb_poly_revert_series(arb_poly_t Qinv, const arb_poly_t Q, slong n, slong prec);

/* Evaluation and interpolation */

void _arb_poly_evaluate_horner(arb_t res, arb_srcptr f, slong len, const arb_t a, slong prec);
void arb_poly_evaluate_horner(arb_t res, const arb_poly_t f, const arb_t a, slong prec);

void _arb_poly_evaluate_rectangular(arb_t y, arb_srcptr poly, slong len, const arb_t x, slong prec);
void arb_poly_evaluate_rectangular(arb_t res, const arb_poly_t f, const arb_t a, slong prec);

void _arb_poly_evaluate(arb_t res, arb_srcptr f, slong len, const arb_t a, slong prec);
void arb_poly_evaluate(arb_t res, const arb_poly_t f, const arb_t a, slong prec);

void _arb_poly_evaluate2_horner(arb_t y, arb_t z, arb_srcptr f, slong len, const arb_t x, slong prec);
void arb_poly_evaluate2_horner(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, slong prec);

void _arb_poly_evaluate2_rectangular(arb_t y, arb_t z, arb_srcptr f, slong len, const arb_t x, slong prec);
void arb_poly_evaluate2_rectangular(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, slong prec);

void _arb_poly_evaluate2(arb_t y, arb_t z, arb_srcptr f, slong len, const arb_t x, slong prec);
void arb_poly_evaluate2(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, slong prec);




void _arb_poly_evaluate_vec_iter(arb_ptr ys, arb_srcptr poly, slong plen,
    arb_srcptr xs, slong n, slong prec);

void arb_poly_evaluate_vec_iter(arb_ptr ys,
        const arb_poly_t poly, arb_srcptr xs, slong n, slong prec);

void _arb_poly_evaluate_vec_fast_precomp(arb_ptr vs, arb_srcptr poly,
    slong plen, arb_ptr * tree, slong len, slong prec);

void _arb_poly_evaluate_vec_fast(arb_ptr ys, arb_srcptr poly, slong plen,
    arb_srcptr xs, slong n, slong prec);

void arb_poly_evaluate_vec_fast(arb_ptr ys,
        const arb_poly_t poly, arb_srcptr xs, slong n, slong prec);

void _arb_poly_interpolate_newton(arb_ptr poly, arb_srcptr xs,
    arb_srcptr ys, slong n, slong prec);

void arb_poly_interpolate_newton(arb_poly_t poly,
    arb_srcptr xs, arb_srcptr ys, slong n, slong prec);

void
_arb_poly_interpolate_barycentric(arb_ptr poly,
    arb_srcptr xs, arb_srcptr ys, slong n, slong prec);

void arb_poly_interpolate_barycentric(arb_poly_t poly,
    arb_srcptr xs, arb_srcptr ys, slong n, slong prec);

void _arb_poly_interpolation_weights(arb_ptr w,
    arb_ptr * tree, slong len, slong prec);

void _arb_poly_interpolate_fast_precomp(arb_ptr poly,
    arb_srcptr ys, arb_ptr * tree, arb_srcptr weights,
    slong len, slong prec);

void _arb_poly_interpolate_fast(arb_ptr poly,
    arb_srcptr xs, arb_srcptr ys, slong len, slong prec);

void arb_poly_interpolate_fast(arb_poly_t poly,
        arb_srcptr xs, arb_srcptr ys, slong n, slong prec);

/* Derivative and integral */

void _arb_poly_derivative(arb_ptr res, arb_srcptr poly, slong len, slong prec);

void arb_poly_derivative(arb_poly_t res, const arb_poly_t poly, slong prec);

void _arb_poly_nth_derivative(arb_ptr res, arb_srcptr poly, ulong n, slong len, slong prec);

void arb_poly_nth_derivative(arb_poly_t res, const arb_poly_t poly, ulong n, slong prec);

void _arb_poly_integral(arb_ptr res, arb_srcptr poly, slong len, slong prec);

void arb_poly_integral(arb_poly_t res, const arb_poly_t poly, slong prec);

/* Transforms */

void arb_poly_borel_transform(arb_poly_t res, const arb_poly_t poly, slong prec);

void _arb_poly_borel_transform(arb_ptr res, arb_srcptr poly, slong len, slong prec);

void arb_poly_inv_borel_transform(arb_poly_t res, const arb_poly_t poly, slong prec);

void _arb_poly_inv_borel_transform(arb_ptr res, arb_srcptr poly, slong len, slong prec);

void _arb_poly_binomial_transform_basecase(arb_ptr b, arb_srcptr a, slong alen, slong len, slong prec);

void arb_poly_binomial_transform_basecase(arb_poly_t b, const arb_poly_t a, slong len, slong prec);

void _arb_poly_binomial_transform_convolution(arb_ptr b, arb_srcptr a, slong alen, slong len, slong prec);

void arb_poly_binomial_transform_convolution(arb_poly_t b, const arb_poly_t a, slong len, slong prec);

void _arb_poly_binomial_transform(arb_ptr b, arb_srcptr a, slong alen, slong len, slong prec);

void arb_poly_binomial_transform(arb_poly_t b, const arb_poly_t a, slong len, slong prec);

void _arb_poly_graeffe_transform(arb_ptr b, arb_srcptr a, slong len, slong prec);

void arb_poly_graeffe_transform(arb_poly_t b, const arb_poly_t a, slong prec);

/* Special functions */

void _arb_poly_pow_ui_trunc_binexp(arb_ptr res,
    arb_srcptr f, slong flen, ulong exp, slong len, slong prec);

void arb_poly_pow_ui_trunc_binexp(arb_poly_t res,
    const arb_poly_t poly, ulong exp, slong len, slong prec);

void _arb_poly_pow_ui(arb_ptr res, arb_srcptr f, slong flen, ulong exp, slong prec);

void arb_poly_pow_ui(arb_poly_t res, const arb_poly_t poly, ulong exp, slong prec);

void _arb_poly_pow_series(arb_ptr h,
    arb_srcptr f, slong flen,
    arb_srcptr g, slong glen, slong len, slong prec);

void arb_poly_pow_series(arb_poly_t h,
    const arb_poly_t f, const arb_poly_t g, slong len, slong prec);

void _arb_poly_pow_arb_series(arb_ptr h,
    arb_srcptr f, slong flen, const arb_t g, slong len, slong prec);

void arb_poly_pow_arb_series(arb_poly_t h,
    const arb_poly_t f, const arb_t g, slong len, slong prec);

void _arb_poly_binomial_pow_arb_series(arb_ptr h, arb_srcptr f, slong flen, const arb_t g, slong len, slong prec);

void _arb_poly_rsqrt_series(arb_ptr g,
    arb_srcptr h, slong hlen, slong len, slong prec);

void arb_poly_rsqrt_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec);

void _arb_poly_sqrt_series(arb_ptr g,
    arb_srcptr h, slong hlen, slong len, slong prec);

void arb_poly_sqrt_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec);

void _arb_poly_log_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec);
void arb_poly_log_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec);

void _arb_poly_log1p_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec);
void arb_poly_log1p_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec);

void _arb_poly_atan_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec);

void arb_poly_atan_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec);

void _arb_poly_asin_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec);

void arb_poly_asin_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec);

void _arb_poly_acos_series(arb_ptr res, arb_srcptr f, slong flen, slong n, slong prec);

void arb_poly_acos_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec);

void _arb_poly_exp_series_basecase(arb_ptr f, arb_srcptr h, slong hlen, slong n, slong prec);
void arb_poly_exp_series_basecase(arb_poly_t f, const arb_poly_t h, slong n, slong prec);
void _arb_poly_exp_series(arb_ptr f, arb_srcptr h, slong hlen, slong n, slong prec);
void arb_poly_exp_series(arb_poly_t f, const arb_poly_t h, slong n, slong prec);

void _arb_poly_sinh_cosh_series_basecase(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec);
void arb_poly_sinh_cosh_series_basecase(arb_poly_t s, arb_poly_t c, const arb_poly_t h, slong n, slong prec);
void _arb_poly_sinh_cosh_series_exponential(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec);
void arb_poly_sinh_cosh_series_exponential(arb_poly_t s, arb_poly_t c, const arb_poly_t h, slong n, slong prec);
void _arb_poly_sinh_cosh_series(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec);
void arb_poly_sinh_cosh_series(arb_poly_t s, arb_poly_t c, const arb_poly_t h, slong n, slong prec);

void _arb_poly_sinh_series(arb_ptr s, arb_srcptr h, slong hlen, slong n, slong prec);
void arb_poly_sinh_series(arb_poly_t s, const arb_poly_t h, slong n, slong prec);

void _arb_poly_cosh_series(arb_ptr c, arb_srcptr h, slong hlen, slong n, slong prec);
void arb_poly_cosh_series(arb_poly_t c, const arb_poly_t h, slong n, slong prec);

void _arb_poly_sin_cos_series(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, slong len, slong prec);
void arb_poly_sin_cos_series(arb_poly_t s, arb_poly_t c, const arb_poly_t h, slong n, slong prec);
void _arb_poly_sin_cos_pi_series(arb_ptr s, arb_ptr c, arb_srcptr h, slong hlen, slong len, slong prec);
void arb_poly_sin_cos_pi_series(arb_poly_t s, arb_poly_t c, const arb_poly_t h, slong n, slong prec);

void _arb_poly_sin_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec);

void arb_poly_sin_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec);

void _arb_poly_cos_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec);

void arb_poly_cos_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec);

void _arb_poly_sin_pi_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec);

void arb_poly_sin_pi_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec);

void _arb_poly_cos_pi_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec);

void arb_poly_cos_pi_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec);

void _arb_poly_cot_pi_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec);

void arb_poly_cot_pi_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec);

void _arb_poly_tan_series(arb_ptr g, arb_srcptr h, slong hlen, slong len, slong prec);

void arb_poly_tan_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec);

void _arb_poly_sinc_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec);
void arb_poly_sinc_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec);

void _arb_poly_sinc_pi_series(arb_ptr g, arb_srcptr h, slong hlen, slong n, slong prec);
void arb_poly_sinc_pi_series(arb_poly_t g, const arb_poly_t h, slong n, slong prec);

void _arb_poly_evaluate_acb_horner(acb_t res, arb_srcptr f, slong len, const acb_t x, slong prec);
void arb_poly_evaluate_acb_horner(acb_t res, const arb_poly_t f, const acb_t a, slong prec);

void _arb_poly_evaluate_acb_rectangular(acb_t y, arb_srcptr poly, slong len, const acb_t x, slong prec);
void arb_poly_evaluate_acb_rectangular(acb_t res, const arb_poly_t f, const acb_t a, slong prec);

void _arb_poly_evaluate_acb(acb_t res, arb_srcptr f, slong len, const acb_t x, slong prec);
void arb_poly_evaluate_acb(acb_t res, const arb_poly_t f, const acb_t a, slong prec);

void _arb_poly_evaluate2_acb_horner(acb_t y, acb_t z, arb_srcptr f, slong len, const acb_t x, slong prec);
void arb_poly_evaluate2_acb_horner(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, slong prec);

void _arb_poly_evaluate2_acb_rectangular(acb_t y, acb_t z, arb_srcptr f, slong len, const acb_t x, slong prec);
void arb_poly_evaluate2_acb_rectangular(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, slong prec);

void _arb_poly_evaluate2_acb(acb_t y, acb_t z, arb_srcptr f, slong len, const acb_t x, slong prec);
void arb_poly_evaluate2_acb(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, slong prec);

void _arb_poly_lambertw_series(arb_ptr res, arb_srcptr z, slong zlen, int flags, slong len, slong prec);
void arb_poly_lambertw_series(arb_poly_t res, const arb_poly_t z, int flags, slong len, slong prec);

void _arb_poly_gamma_series(arb_ptr res, arb_srcptr h, slong hlen, slong len, slong prec);
void arb_poly_gamma_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec);

void _arb_poly_rgamma_series(arb_ptr res, arb_srcptr h, slong hlen, slong len, slong prec);
void arb_poly_rgamma_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec);

void _arb_poly_lgamma_series(arb_ptr res, arb_srcptr h, slong hlen, slong len, slong prec);
void arb_poly_lgamma_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec);

void _arb_poly_digamma_series(arb_ptr res, arb_srcptr h, slong hlen, slong len, slong prec);
void arb_poly_digamma_series(arb_poly_t res, const arb_poly_t f, slong n, slong prec);

void _arb_poly_rising_ui_series(arb_ptr res, arb_srcptr f, slong flen, ulong r, slong trunc, slong prec);
void arb_poly_rising_ui_series(arb_poly_t res, const arb_poly_t f, ulong r, slong trunc, slong prec);

void _arb_poly_zeta_series(arb_ptr res, arb_srcptr h, slong hlen, const arb_t a, int deflate, slong len, slong prec);
void arb_poly_zeta_series(arb_poly_t res, const arb_poly_t f, const arb_t a, int deflate, slong n, slong prec);

void _arb_poly_riemann_siegel_theta_series(arb_ptr res, arb_srcptr h, slong hlen, slong len, slong prec);
void arb_poly_riemann_siegel_theta_series(arb_poly_t res, const arb_poly_t h, slong n, slong prec);

void _arb_poly_riemann_siegel_z_series(arb_ptr res, arb_srcptr h, slong hlen, slong len, slong prec);
void arb_poly_riemann_siegel_z_series(arb_poly_t res, const arb_poly_t h, slong n, slong prec);

slong _arb_poly_swinnerton_dyer_ui_prec(ulong n);
void _arb_poly_swinnerton_dyer_ui(arb_ptr T, ulong n, slong trunc, slong prec);
void arb_poly_swinnerton_dyer_ui(arb_poly_t poly, ulong n, slong prec);

/* Root-finding */

void _arb_poly_newton_convergence_factor(arf_t convergence_factor,
    arb_srcptr poly, slong len,
    const arb_t convergence_interval, slong prec);

int _arb_poly_newton_step(arb_t xnew, arb_srcptr poly, slong len,
    const arb_t x,
    const arb_t convergence_interval,
    const arf_t convergence_factor, slong prec);

void _arb_poly_newton_refine_root(arb_t r, arb_srcptr poly,
    slong len, const arb_t start,
    const arb_t convergence_interval,
    const arf_t convergence_factor,
    slong eval_extra_prec,
    slong prec);

void _arb_poly_root_bound_fujiwara(mag_t bound, arb_srcptr poly, slong len);

void arb_poly_root_bound_fujiwara(mag_t bound, arb_poly_t poly);

ARB_POLY_INLINE slong
arb_poly_allocated_bytes(const arb_poly_t x)
{
    return _arb_vec_allocated_bytes(x->coeffs, x->alloc);
}

/* Macros */


/* counts zero bits in the binary representation of e */
ARB_POLY_INLINE int
n_zerobits(mp_limb_t e)
{
    int zeros = 0;

    while (e > 1)
    {
        zeros += !(e & 1);
        e >>= 1;
    }

    return zeros;
}

/* Computes the length of the result when raising a polynomial of
   length *len* to the power *exp* and truncating to length *trunc*,
   without overflow. Assumes poly_len >= 1. */
ARB_POLY_INLINE slong
poly_pow_length(slong poly_len, ulong exp, slong trunc)
{
    mp_limb_t hi, lo;
    umul_ppmm(hi, lo, poly_len - 1, exp);
    add_ssaaaa(hi, lo, hi, lo, 0, 1);
    if (hi != 0 || lo > (mp_limb_t) WORD_MAX)
        return trunc;
    return FLINT_MIN((slong) lo, trunc);
}

#ifndef NEWTON_INIT

#define NEWTON_INIT(from, to) \
    { \
        slong __steps[FLINT_BITS], __i, __from, __to; \
        __steps[__i = 0] = __to = (to); \
        __from = (from); \
        while (__to > __from) \
            __steps[++__i] = (__to = (__to + 1) / 2); \

#define NEWTON_BASECASE(bc_to) { slong bc_to = __to;

#define NEWTON_END_BASECASE }

#define NEWTON_LOOP(step_from, step_to) \
        { \
            for (__i--; __i >= 0; __i--) \
            { \
                slong step_from = __steps[__i+1]; \
                slong step_to = __steps[__i]; \

#define NEWTON_END_LOOP }}

#define NEWTON_END }

#endif

#ifdef __cplusplus
}
#endif

#endif

