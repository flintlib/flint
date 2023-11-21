/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_POLY_H
#define ACB_POLY_H

#ifdef ACB_POLY_INLINES_C
#define ACB_POLY_INLINE
#else
#define ACB_POLY_INLINE static inline
#endif

#include "fmpq_types.h"
#include "acb.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Memory management */

void acb_poly_init(acb_poly_t poly);

void acb_poly_init2(acb_poly_t poly, slong len);

void acb_poly_clear(acb_poly_t poly);

void acb_poly_fit_length(acb_poly_t poly, slong len);

void _acb_poly_set_length(acb_poly_t poly, slong len);

void _acb_poly_normalise(acb_poly_t poly);

ACB_POLY_INLINE void
acb_poly_swap(acb_poly_t poly1, acb_poly_t poly2)
{
    acb_poly_struct t = *poly1;
    *poly1 = *poly2;
    *poly2 = t;
}

ACB_POLY_INLINE slong acb_poly_length(const acb_poly_t poly)
{
    return poly->length;
}

ACB_POLY_INLINE slong acb_poly_degree(const acb_poly_t poly)
{
    return poly->length - 1;
}

slong acb_poly_valuation(const acb_poly_t poly);

ACB_POLY_INLINE int
acb_poly_is_zero(const acb_poly_t z)
{
    return acb_poly_length(z) == 0;
}

ACB_POLY_INLINE int
acb_poly_is_one(const acb_poly_t z)
{
    return (acb_poly_length(z) == 1) && acb_is_one(z->coeffs);
}

ACB_POLY_INLINE int
acb_poly_is_x(const acb_poly_t z)
{
    return (acb_poly_length(z) == 2) && acb_is_zero(z->coeffs)
        && acb_is_one(z->coeffs + 1);
}

ACB_POLY_INLINE void acb_poly_zero(acb_poly_t poly)
{
    poly->length = 0;
}

ACB_POLY_INLINE void
acb_poly_one(acb_poly_t poly)
{
    acb_poly_fit_length(poly, 1);
    acb_one(poly->coeffs);
    _acb_poly_set_length(poly, 1);
}

void acb_poly_set_coeff_si(acb_poly_t poly, slong n, slong x);

void acb_poly_set_coeff_acb(acb_poly_t poly, slong n, const acb_t x);

void acb_poly_get_coeff_acb(acb_t x, const acb_poly_t poly, slong n);

#define acb_poly_get_coeff_ptr(poly, n) \
    ((n) < (poly)->length ? (poly)->coeffs + (n) : NULL)

void _acb_poly_shift_right(acb_ptr res, acb_srcptr poly, slong len, slong n);

void acb_poly_shift_right(acb_poly_t res, const acb_poly_t poly, slong n);

void _acb_poly_shift_left(acb_ptr res, acb_srcptr poly, slong len, slong n);

void acb_poly_shift_left(acb_poly_t res, const acb_poly_t poly, slong n);

ACB_POLY_INLINE void
acb_poly_truncate(acb_poly_t poly, slong newlen)
{
    if (poly->length > newlen)
    {
        slong i;
        for (i = newlen; i < poly->length; i++)
            acb_zero(poly->coeffs + i);
        poly->length = newlen;
        _acb_poly_normalise(poly);
    }
}

void _acb_poly_majorant(arb_ptr res, acb_srcptr vec, slong len, slong prec);

void acb_poly_majorant(arb_poly_t res, const acb_poly_t poly, slong prec);

#ifdef FLINT_HAVE_FILE
void acb_poly_fprintd(FILE * file, const acb_poly_t poly, slong digits);
#endif

void acb_poly_printd(const acb_poly_t poly, slong digits);

void _acb_poly_evaluate_horner(acb_t res, acb_srcptr f, slong len, const acb_t a, slong prec);
void acb_poly_evaluate_horner(acb_t res, const acb_poly_t f, const acb_t a, slong prec);

void _acb_poly_evaluate_rectangular(acb_t y, acb_srcptr poly, slong len, const acb_t x, slong prec);
void acb_poly_evaluate_rectangular(acb_t res, const acb_poly_t f, const acb_t a, slong prec);

void _acb_poly_evaluate(acb_t res, acb_srcptr f, slong len, const acb_t a, slong prec);
void acb_poly_evaluate(acb_t res, const acb_poly_t f, const acb_t a, slong prec);

void _acb_poly_evaluate2_horner(acb_t y, acb_t z, acb_srcptr f, slong len, const acb_t x, slong prec);
void acb_poly_evaluate2_horner(acb_t y, acb_t z, const acb_poly_t f, const acb_t x, slong prec);

void _acb_poly_evaluate2_rectangular(acb_t y, acb_t z, acb_srcptr f, slong len, const acb_t x, slong prec);
void acb_poly_evaluate2_rectangular(acb_t y, acb_t z, const acb_poly_t f, const acb_t x, slong prec);

void _acb_poly_evaluate2(acb_t y, acb_t z, acb_srcptr f, slong len, const acb_t x, slong prec);
void acb_poly_evaluate2(acb_t y, acb_t z, const acb_poly_t f, const acb_t x, slong prec);

void _acb_poly_derivative(acb_ptr res, acb_srcptr poly, slong len, slong prec);

void acb_poly_derivative(acb_poly_t res, const acb_poly_t poly, slong prec);

void _acb_poly_nth_derivative(acb_ptr res, acb_srcptr poly, ulong n, slong len, slong prec);

void acb_poly_nth_derivative(acb_poly_t res, const acb_poly_t poly, ulong n, slong prec);

void _acb_poly_integral(acb_ptr res, acb_srcptr poly, slong len, slong prec);

void acb_poly_integral(acb_poly_t res, const acb_poly_t poly, slong prec);

/* Transforms */

void acb_poly_borel_transform(acb_poly_t res, const acb_poly_t poly, slong prec);

void _acb_poly_borel_transform(acb_ptr res, acb_srcptr poly, slong len, slong prec);

void acb_poly_inv_borel_transform(acb_poly_t res, const acb_poly_t poly, slong prec);

void _acb_poly_inv_borel_transform(acb_ptr res, acb_srcptr poly, slong len, slong prec);

void _acb_poly_binomial_transform_basecase(acb_ptr b, acb_srcptr a, slong alen, slong len, slong prec);

void acb_poly_binomial_transform_basecase(acb_poly_t b, const acb_poly_t a, slong len, slong prec);

void _acb_poly_binomial_transform_convolution(acb_ptr b, acb_srcptr a, slong alen, slong len, slong prec);

void acb_poly_binomial_transform_convolution(acb_poly_t b, const acb_poly_t a, slong len, slong prec);

void _acb_poly_binomial_transform(acb_ptr b, acb_srcptr a, slong alen, slong len, slong prec);

void acb_poly_binomial_transform(acb_poly_t b, const acb_poly_t a, slong len, slong prec);

void _acb_poly_graeffe_transform(acb_ptr b, acb_srcptr a, slong len, slong prec);

void acb_poly_graeffe_transform(acb_poly_t b, const acb_poly_t a, slong prec);



void acb_poly_set(acb_poly_t dest, const acb_poly_t src);

void acb_poly_set_round(acb_poly_t dest, const acb_poly_t src, slong prec);

void acb_poly_set_trunc(acb_poly_t res, const acb_poly_t poly, slong n);

void acb_poly_set_trunc_round(acb_poly_t res, const acb_poly_t poly, slong n, slong prec);

void acb_poly_set_arb_poly(acb_poly_t poly, const arb_poly_t re);

void acb_poly_set2_arb_poly(acb_poly_t poly, const arb_poly_t re, const arb_poly_t im);

void acb_poly_set_fmpq_poly(acb_poly_t poly, const fmpq_poly_t re, slong prec);

void acb_poly_set2_fmpq_poly(acb_poly_t poly, const fmpq_poly_t re, const fmpq_poly_t im, slong prec);

void acb_poly_set_fmpz_poly(acb_poly_t poly, const fmpz_poly_t src, slong prec);

void acb_poly_set2_fmpz_poly(acb_poly_t poly, const fmpz_poly_t re, const fmpz_poly_t im, slong prec);

int acb_poly_get_unique_fmpz_poly(fmpz_poly_t res, const acb_poly_t src);

ACB_POLY_INLINE void
acb_poly_set_acb(acb_poly_t poly, const acb_t c)
{
    acb_poly_fit_length(poly, 1);
    acb_set(poly->coeffs, c);
    _acb_poly_set_length(poly, !acb_is_zero(poly->coeffs));
}

void acb_poly_set_si(acb_poly_t poly, slong c);

void acb_poly_randtest(acb_poly_t poly, flint_rand_t state, slong len, slong prec, slong mag_bits);

int acb_poly_equal(const acb_poly_t A, const acb_poly_t B);

int acb_poly_contains_fmpz_poly(const acb_poly_t poly1, const fmpz_poly_t poly2);

int acb_poly_contains_fmpq_poly(const acb_poly_t poly1, const fmpq_poly_t poly2);

int _acb_poly_overlaps(acb_srcptr poly1, slong len1,
        acb_srcptr poly2, slong len2);

int acb_poly_overlaps(const acb_poly_t poly1, const acb_poly_t poly2);

int acb_poly_contains(const acb_poly_t poly1, const acb_poly_t poly2);

ACB_POLY_INLINE int
acb_poly_is_real(const acb_poly_t poly)
{
    return _acb_vec_is_real(poly->coeffs, poly->length);
}

void _acb_poly_add(acb_ptr res, acb_srcptr poly1, slong len1,
    acb_srcptr poly2, slong len2, slong prec);

void acb_poly_add(acb_poly_t res, const acb_poly_t poly1,
              const acb_poly_t poly2, slong prec);

void acb_poly_add_si(acb_poly_t res, const acb_poly_t poly, slong c, slong prec);

void _acb_poly_sub(acb_ptr res, acb_srcptr poly1, slong len1,
    acb_srcptr poly2, slong len2, slong prec);

void acb_poly_sub(acb_poly_t res, const acb_poly_t poly1,
              const acb_poly_t poly2, slong prec);

void acb_poly_add_series(acb_poly_t res, const acb_poly_t poly1,
              const acb_poly_t poly2, slong len, slong prec);

void acb_poly_sub_series(acb_poly_t res, const acb_poly_t poly1,
              const acb_poly_t poly2, slong len, slong prec);

ACB_POLY_INLINE void
acb_poly_neg(acb_poly_t res, const acb_poly_t poly)
{
    acb_poly_fit_length(res, poly->length);
    _acb_vec_neg(res->coeffs, poly->coeffs, poly->length);
    _acb_poly_set_length(res, poly->length);
}

ACB_POLY_INLINE void
acb_poly_scalar_mul_2exp_si(acb_poly_t res, const acb_poly_t poly, slong c)
{
    acb_poly_fit_length(res, poly->length);
    _acb_vec_scalar_mul_2exp_si(res->coeffs, poly->coeffs, poly->length, c);
    _acb_poly_set_length(res, poly->length);
}

ACB_POLY_INLINE void
acb_poly_scalar_mul(acb_poly_t res, const acb_poly_t poly, const acb_t c, slong prec)
{
    acb_poly_fit_length(res, poly->length);
    _acb_vec_scalar_mul(res->coeffs, poly->coeffs, poly->length, c, prec);
    _acb_poly_set_length(res, poly->length);
    _acb_poly_normalise(res);
}

ACB_POLY_INLINE void
acb_poly_scalar_div(acb_poly_t res, const acb_poly_t poly, const acb_t c, slong prec)
{
    acb_poly_fit_length(res, poly->length);
    _acb_vec_scalar_div(res->coeffs, poly->coeffs, poly->length, c, prec);
    _acb_poly_set_length(res, poly->length);
    _acb_poly_normalise(res);
}

void acb_poly_mullow_classical(acb_poly_t res, const acb_poly_t poly1,
                                            const acb_poly_t poly2,
                                                slong n, slong prec);

void _acb_poly_mullow_classical(acb_ptr res,
    acb_srcptr poly1, slong len1,
    acb_srcptr poly2, slong len2, slong n, slong prec);

void _acb_poly_mullow_transpose(acb_ptr res,
    acb_srcptr poly1, slong len1,
    acb_srcptr poly2, slong len2, slong n, slong prec);

void acb_poly_mullow_transpose(acb_poly_t res, const acb_poly_t poly1,
                                            const acb_poly_t poly2,
                                                slong n, slong prec);

void _acb_poly_mullow_transpose_gauss(acb_ptr res,
    acb_srcptr poly1, slong len1,
    acb_srcptr poly2, slong len2, slong n, slong prec);

void acb_poly_mullow_transpose_gauss(acb_poly_t res, const acb_poly_t poly1,
                                            const acb_poly_t poly2,
                                                slong n, slong prec);

void _acb_poly_mullow(acb_ptr res,
    acb_srcptr poly1, slong len1,
    acb_srcptr poly2, slong len2, slong n, slong prec);

void acb_poly_mullow(acb_poly_t res, const acb_poly_t poly1,
                                            const acb_poly_t poly2,
                                                slong n, slong prec);

void _acb_poly_mul(acb_ptr C,
    acb_srcptr A, slong lenA,
    acb_srcptr B, slong lenB, slong prec);

void acb_poly_mul(acb_poly_t res, const acb_poly_t poly1,
              const acb_poly_t poly2, slong prec);

ACB_POLY_INLINE void
_acb_poly_mul_monic(acb_ptr res, acb_srcptr poly1, slong len1,
    acb_srcptr poly2, slong len2, slong prec)
{
    if (len1 + len2 - 2 > 0)
        _acb_poly_mullow(res, poly1, len1, poly2, len2, len1 + len2 - 2, prec);
    acb_one(res + len1 + len2 - 2);
}

void _acb_poly_inv_series(acb_ptr Qinv, acb_srcptr Q, slong Qlen, slong len, slong prec);

void acb_poly_inv_series(acb_poly_t Qinv, const acb_poly_t Q, slong n, slong prec);

void  _acb_poly_div_series(acb_ptr Q, acb_srcptr A, slong Alen,
    acb_srcptr B, slong Blen, slong n, slong prec);

void acb_poly_div_series(acb_poly_t Q, const acb_poly_t A, const acb_poly_t B, slong n, slong prec);

void _acb_poly_reverse(acb_ptr res, acb_srcptr poly, slong len, slong n);

void _acb_poly_div(acb_ptr Q,
    acb_srcptr A, slong lenA,
    acb_srcptr B, slong lenB, slong prec);

void _acb_poly_divrem(acb_ptr Q, acb_ptr R,
    acb_srcptr A, slong lenA,
    acb_srcptr B, slong lenB, slong prec);

void _acb_poly_rem(acb_ptr R,
    acb_srcptr A, slong lenA,
    acb_srcptr B, slong lenB, slong prec);

int acb_poly_divrem(acb_poly_t Q, acb_poly_t R,
                             const acb_poly_t A, const acb_poly_t B, slong prec);

void _acb_poly_div_root(acb_ptr Q, acb_t R, acb_srcptr A,
    slong len, const acb_t c, slong prec);

/* Composition */

void _acb_poly_taylor_shift(acb_ptr poly, const acb_t c, slong n, slong prec);
void acb_poly_taylor_shift(acb_poly_t g, const acb_poly_t f, const acb_t c, slong prec);
void _acb_poly_compose(acb_ptr res, acb_srcptr poly1, slong len1, acb_srcptr poly2, slong len2, slong prec);
void acb_poly_compose(acb_poly_t res, const acb_poly_t poly1, const acb_poly_t poly2, slong prec);

void _acb_poly_compose_series(acb_ptr res, acb_srcptr poly1, slong len1,
                            acb_srcptr poly2, slong len2, slong n, slong prec);

void acb_poly_compose_series(acb_poly_t res,
                    const acb_poly_t poly1,
                    const acb_poly_t poly2, slong n, slong prec);

/* Reversion */

void _acb_poly_revert_series(acb_ptr Qinv, acb_srcptr Q, slong Qlen, slong n, slong prec);
void acb_poly_revert_series(acb_poly_t Qinv, const acb_poly_t Q, slong n, slong prec);



void
_acb_poly_evaluate_vec_fast_precomp(acb_ptr vs, acb_srcptr poly,
    slong plen, acb_ptr * tree, slong len, slong prec);

void _acb_poly_evaluate_vec_fast(acb_ptr ys, acb_srcptr poly, slong plen,
    acb_srcptr xs, slong n, slong prec);

void
acb_poly_evaluate_vec_fast(acb_ptr ys,
        const acb_poly_t poly, acb_srcptr xs, slong n, slong prec);

void
_acb_poly_evaluate_vec_iter(acb_ptr ys, acb_srcptr poly, slong plen,
    acb_srcptr xs, slong n, slong prec);

void
acb_poly_evaluate_vec_iter(acb_ptr ys,
        const acb_poly_t poly, acb_srcptr xs, slong n, slong prec);

void
_acb_poly_interpolate_barycentric(acb_ptr poly,
    acb_srcptr xs, acb_srcptr ys, slong n, slong prec);

void
acb_poly_interpolate_barycentric(acb_poly_t poly,
    acb_srcptr xs, acb_srcptr ys, slong n, slong prec);

void
_acb_poly_interpolation_weights(acb_ptr w,
    acb_ptr * tree, slong len, slong prec);

void
_acb_poly_interpolate_fast_precomp(acb_ptr poly,
    acb_srcptr ys, acb_ptr * tree, acb_srcptr weights,
    slong len, slong prec);

void
_acb_poly_interpolate_fast(acb_ptr poly,
    acb_srcptr xs, acb_srcptr ys, slong len, slong prec);

void
acb_poly_interpolate_fast(acb_poly_t poly,
        acb_srcptr xs, acb_srcptr ys, slong n, slong prec);

void
_acb_poly_interpolate_newton(acb_ptr poly, acb_srcptr xs,
    acb_srcptr ys, slong n, slong prec);

void
acb_poly_interpolate_newton(acb_poly_t poly,
    acb_srcptr xs, acb_srcptr ys, slong n, slong prec);

void
_acb_poly_product_roots(acb_ptr poly, acb_srcptr xs, slong n, slong prec);

void
acb_poly_product_roots(acb_poly_t poly, acb_srcptr xs, slong n, slong prec);

acb_ptr * _acb_poly_tree_alloc(slong len);

void _acb_poly_tree_free(acb_ptr * tree, slong len);

void
_acb_poly_tree_build(acb_ptr * tree, acb_srcptr roots, slong len, slong prec);


void _acb_poly_root_inclusion(acb_t r, const acb_t m,
    acb_srcptr poly,
    acb_srcptr polyder, slong len, slong prec);

slong _acb_poly_validate_roots(acb_ptr roots,
        acb_srcptr poly, slong len, slong prec);

void _acb_poly_refine_roots_durand_kerner(acb_ptr roots,
        acb_srcptr poly, slong len, slong prec);

slong _acb_poly_find_roots(acb_ptr roots,
    acb_srcptr poly,
    acb_srcptr initial, slong len, slong maxiter, slong prec);

slong acb_poly_find_roots(acb_ptr roots,
    const acb_poly_t poly, acb_srcptr initial,
    slong maxiter, slong prec);

void _acb_poly_root_bound_fujiwara(mag_t bound, acb_srcptr poly, slong len);

void acb_poly_root_bound_fujiwara(mag_t bound, acb_poly_t poly);

int _acb_poly_validate_real_roots(acb_srcptr roots, acb_srcptr poly, slong len, slong prec);

int acb_poly_validate_real_roots(acb_srcptr roots, const acb_poly_t poly, slong prec);

/* Special functions */

void _acb_poly_pow_ui_trunc_binexp(acb_ptr res,
    acb_srcptr f, slong flen, ulong exp, slong len, slong prec);

void acb_poly_pow_ui_trunc_binexp(acb_poly_t res,
    const acb_poly_t poly, ulong exp, slong len, slong prec);

void _acb_poly_pow_ui(acb_ptr res, acb_srcptr f, slong flen, ulong exp, slong prec);

void acb_poly_pow_ui(acb_poly_t res, const acb_poly_t poly, ulong exp, slong prec);

void _acb_poly_rsqrt_series(acb_ptr g, acb_srcptr h, slong hlen, slong len, slong prec);

void acb_poly_rsqrt_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec);

void _acb_poly_sqrt_series(acb_ptr g, acb_srcptr h, slong hlen, slong len, slong prec);

void acb_poly_sqrt_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec);

void _acb_poly_log_series(acb_ptr res, acb_srcptr f, slong flen, slong n, slong prec);
void acb_poly_log_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec);

void _acb_poly_log1p_series(acb_ptr res, acb_srcptr f, slong flen, slong n, slong prec);
void acb_poly_log1p_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec);

void _acb_poly_atan_series(acb_ptr res, acb_srcptr f, slong flen, slong n, slong prec);

void acb_poly_atan_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec);

void _acb_poly_exp_series_basecase(acb_ptr f, acb_srcptr h, slong hlen, slong n, slong prec);
void acb_poly_exp_series_basecase(acb_poly_t f, const acb_poly_t h, slong n, slong prec);
void _acb_poly_exp_series(acb_ptr f, acb_srcptr h, slong hlen, slong n, slong prec);
void acb_poly_exp_series(acb_poly_t f, const acb_poly_t h, slong n, slong prec);

void _acb_poly_exp_pi_i_series(acb_ptr f, acb_srcptr h, slong hlen, slong n, slong prec);
void acb_poly_exp_pi_i_series(acb_poly_t f, const acb_poly_t h, slong n, slong prec);

void _acb_poly_sinh_cosh_series_basecase(acb_ptr s, acb_ptr c, acb_srcptr h, slong hlen, slong n, slong prec);
void acb_poly_sinh_cosh_series_basecase(acb_poly_t s, acb_poly_t c, const acb_poly_t h, slong n, slong prec);
void _acb_poly_sinh_cosh_series_exponential(acb_ptr s, acb_ptr c, acb_srcptr h, slong hlen, slong n, slong prec);
void acb_poly_sinh_cosh_series_exponential(acb_poly_t s, acb_poly_t c, const acb_poly_t h, slong n, slong prec);
void _acb_poly_sinh_cosh_series(acb_ptr s, acb_ptr c, acb_srcptr h, slong hlen, slong n, slong prec);
void acb_poly_sinh_cosh_series(acb_poly_t s, acb_poly_t c, const acb_poly_t h, slong n, slong prec);

void _acb_poly_sinh_series(acb_ptr s, acb_srcptr h, slong hlen, slong n, slong prec);
void acb_poly_sinh_series(acb_poly_t s, const acb_poly_t h, slong n, slong prec);

void _acb_poly_cosh_series(acb_ptr c, acb_srcptr h, slong hlen, slong n, slong prec);
void acb_poly_cosh_series(acb_poly_t c, const acb_poly_t h, slong n, slong prec);

void _acb_poly_sin_cos_series(acb_ptr s, acb_ptr c, acb_srcptr h, slong hlen, slong len, slong prec);
void acb_poly_sin_cos_series(acb_poly_t s, acb_poly_t c, const acb_poly_t h, slong n, slong prec);

void _acb_poly_sin_series(acb_ptr g, acb_srcptr h, slong hlen, slong n, slong prec);

void acb_poly_sin_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec);

void _acb_poly_cos_series(acb_ptr g, acb_srcptr h, slong hlen, slong n, slong prec);

void acb_poly_cos_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec);

void _acb_poly_sin_cos_pi_series(acb_ptr s, acb_ptr c, acb_srcptr h, slong hlen, slong len, slong prec);
void acb_poly_sin_cos_pi_series(acb_poly_t s, acb_poly_t c, const acb_poly_t h, slong n, slong prec);

void _acb_poly_sin_pi_series(acb_ptr g, acb_srcptr h, slong hlen, slong n, slong prec);

void acb_poly_sin_pi_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec);

void _acb_poly_cos_pi_series(acb_ptr g, acb_srcptr h, slong hlen, slong n, slong prec);

void acb_poly_cos_pi_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec);

void _acb_poly_cot_pi_series(acb_ptr g, acb_srcptr h, slong hlen, slong len, slong prec);

void acb_poly_cot_pi_series(acb_poly_t res, const acb_poly_t f, slong len, slong prec);

void _acb_poly_tan_series(acb_ptr g, acb_srcptr h, slong hlen, slong len, slong prec);

void acb_poly_tan_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec);

void _acb_poly_sinc_series(acb_ptr g, acb_srcptr h, slong hlen, slong n, slong prec);
void acb_poly_sinc_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec);

void _acb_poly_lambertw_series(acb_ptr res, acb_srcptr z, slong zlen, const fmpz_t k, int flags, slong len, slong prec);
void acb_poly_lambertw_series(acb_poly_t res, const acb_poly_t z, const fmpz_t k, int flags, slong len, slong prec);

void _acb_poly_gamma_series(acb_ptr res, acb_srcptr h, slong hlen, slong len, slong prec);

void acb_poly_gamma_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec);

void _acb_poly_rgamma_series(acb_ptr res, acb_srcptr h, slong hlen, slong len, slong prec);

void acb_poly_rgamma_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec);

void _acb_poly_lgamma_series(acb_ptr res, acb_srcptr h, slong hlen, slong len, slong prec);

void acb_poly_lgamma_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec);

void _acb_poly_digamma_series(acb_ptr res, acb_srcptr h, slong hlen, slong len, slong prec);

void acb_poly_digamma_series(acb_poly_t res, const acb_poly_t f, slong n, slong prec);

void _acb_poly_rising_ui_series(acb_ptr res, acb_srcptr f, slong flen, ulong r, slong trunc, slong prec);

void acb_poly_rising_ui_series(acb_poly_t res, const acb_poly_t f, ulong r, slong trunc, slong prec);

void _acb_poly_pow_acb_series(acb_ptr h,
    acb_srcptr f, slong flen, const acb_t g, slong len, slong prec);

void acb_poly_pow_acb_series(acb_poly_t h,
    const acb_poly_t f, const acb_t g, slong len, slong prec);

void _acb_poly_pow_series(acb_ptr h,
    acb_srcptr f, slong flen,
    acb_srcptr g, slong glen, slong len, slong prec);

void acb_poly_pow_series(acb_poly_t h,
    const acb_poly_t f, const acb_poly_t g, slong len, slong prec);

void
_acb_poly_binomial_pow_acb_series(acb_ptr h, acb_srcptr f, slong flen,
    const acb_t g, slong len, slong prec);

/* TODO: document */
ACB_POLY_INLINE void
_acb_poly_acb_pow_cpx(acb_ptr w, const acb_t a, const acb_t b, slong len, slong prec)
{
    if (len == 1)
    {
        acb_pow(w, a, b, prec);
    }
    else
    {
        acb_t log_a;
        slong k;

        acb_init(log_a);

        acb_log(log_a, a, prec);
        acb_mul(w, log_a, b, prec);
        acb_exp(w, w, prec);

        for (k = 1; k < len; k++)
        {
            acb_mul(w + k, w + k - 1, log_a, prec);
            acb_div_ui(w + k, w + k, k, prec);
        }

        acb_clear(log_a);
    }
}

#define _acb_poly_pow_cpx _acb_poly_acb_pow_cpx

/* TODO: document */
void _acb_poly_acb_invpow_cpx(acb_ptr res, const acb_t N, const acb_t c, slong trunc, slong prec);
/* TODO: document */
void _acb_poly_mullow_cpx(acb_ptr res, acb_srcptr src, slong len, const acb_t c, slong trunc, slong prec);

void _acb_poly_powsum_series_naive(acb_ptr z, const acb_t s, const acb_t a, const acb_t q, slong n, slong len, slong prec);
void _acb_poly_powsum_series_naive_threaded(acb_ptr z, const acb_t s, const acb_t a, const acb_t q, slong n, slong len, slong prec);
void _acb_poly_powsum_one_series_sieved(acb_ptr z, const acb_t s, slong n, slong len, slong prec);

void _acb_poly_zeta_em_sum(acb_ptr z, const acb_t s, const acb_t a, int deflate, ulong N, ulong M, slong d, slong prec);
void _acb_poly_zeta_em_choose_param(mag_t bound, ulong * N, ulong * M, const acb_t s, const acb_t a, slong d, slong target, slong prec);
void _acb_poly_zeta_em_bound1(mag_t bound, const acb_t s, const acb_t a, slong N, slong M, slong d, slong wp);
void _acb_poly_zeta_em_bound(arb_ptr vec, const acb_t s, const acb_t a, ulong N, ulong M, slong d, slong wp);

void _acb_poly_zeta_em_tail_naive(acb_ptr sum, const acb_t s, const acb_t Na, acb_srcptr Nasx, slong M, slong len, slong prec);
void _acb_poly_zeta_em_tail_bsplit(acb_ptr z, const acb_t s, const acb_t Na, acb_srcptr Nasx, slong M, slong len, slong prec);

void _acb_poly_zeta_cpx_series(acb_ptr z, const acb_t s, const acb_t a, int deflate, slong d, slong prec);

void _acb_poly_zeta_series(acb_ptr res, acb_srcptr h, slong hlen, const acb_t a, int deflate, slong len, slong prec);

void acb_poly_zeta_series(acb_poly_t res, const acb_poly_t f, const acb_t a, int deflate, slong n, slong prec);

void _acb_poly_polylog_cpx_zeta(acb_ptr w, const acb_t s, const acb_t z, slong len, slong prec);
void _acb_poly_polylog_cpx_small(acb_ptr w, const acb_t s, const acb_t z, slong len, slong prec);
void _acb_poly_polylog_cpx(acb_ptr w, const acb_t s, const acb_t z, slong len, slong prec);

void _acb_poly_polylog_series(acb_ptr res, acb_srcptr s, slong slen, const acb_t z, slong len, slong prec);
void acb_poly_polylog_series(acb_poly_t res, const acb_poly_t s, const acb_t z, slong n, slong prec);

void _acb_poly_agm1_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec);
void acb_poly_agm1_series(acb_poly_t res, const acb_poly_t z, slong n, slong prec);
void _acb_poly_elliptic_k_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec);
void acb_poly_elliptic_k_series(acb_poly_t res, const acb_poly_t z, slong n, slong prec);
void _acb_poly_elliptic_p_series(acb_ptr res, acb_srcptr z, slong zlen, const acb_t tau, slong len, slong prec);
void acb_poly_elliptic_p_series(acb_poly_t res, const acb_poly_t z, const acb_t tau, slong n, slong prec);

void _acb_poly_erf_series(acb_ptr g, acb_srcptr h, slong hlen, slong n, slong prec);
void acb_poly_erf_series(acb_poly_t g, const acb_poly_t h, slong n, slong prec);

ACB_POLY_INLINE slong
acb_poly_allocated_bytes(const acb_poly_t x)
{
    return _acb_vec_allocated_bytes(x->coeffs, x->alloc);
}

#ifdef __cplusplus
}
#endif

#endif
