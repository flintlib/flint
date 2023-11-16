/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ARB_H
#define ARB_H

#ifdef ARB_INLINES_C
#define ARB_INLINE
#else
#define ARB_INLINE static inline
#endif

#include "arf.h"
#include "arb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define arb_midref(x) (&(x)->mid)
#define arb_radref(x) (&(x)->rad)

#define ARB_IS_LAGOM(x) (ARF_IS_LAGOM(arb_midref(x)) && MAG_IS_LAGOM(arb_radref(x)))

#define ARB_RND ARF_RND_DOWN

ARB_INLINE void
arb_init(arb_t x)
{
    arf_init(arb_midref(x));
    mag_init(arb_radref(x));
}

void arb_clear(arb_t x);

arb_ptr _arb_vec_init(slong n);
void _arb_vec_clear(arb_ptr v, slong n);

ARB_INLINE arf_ptr arb_mid_ptr(arb_t z) { return arb_midref(z); }
ARB_INLINE mag_ptr arb_rad_ptr(arb_t z) { return arb_radref(z); }

ARB_INLINE int
arb_is_exact(const arb_t x)
{
    return mag_is_zero(arb_radref(x));
}

ARB_INLINE int
arb_equal(const arb_t x, const arb_t y)
{
    return arf_equal(arb_midref(x), arb_midref(y)) &&
           mag_equal(arb_radref(x), arb_radref(y));
}

ARB_INLINE int
arb_equal_si(const arb_t x, slong y)
{
    return arf_equal_si(arb_midref(x), y) && mag_is_zero(arb_radref(x));
}

/* implementations are in arb/richcmp.c */
int arb_eq(const arb_t x, const arb_t y);
int arb_ne(const arb_t x, const arb_t y);
int arb_lt(const arb_t x, const arb_t y);
int arb_le(const arb_t x, const arb_t y);
int arb_gt(const arb_t x, const arb_t y);
int arb_ge(const arb_t x, const arb_t y);

void arb_zero(arb_t x);

ARB_INLINE int
arb_is_zero(const arb_t x)
{
    return arf_is_zero(arb_midref(x)) && mag_is_zero(arb_radref(x));
}

ARB_INLINE void
arb_pos_inf(arb_t x)
{
    arf_pos_inf(arb_midref(x));
    mag_zero(arb_radref(x));
}

ARB_INLINE void
arb_neg_inf(arb_t x)
{
    arf_neg_inf(arb_midref(x));
    mag_zero(arb_radref(x));
}

ARB_INLINE void
arb_zero_pm_inf(arb_t x)
{
    arf_zero(arb_midref(x));
    mag_inf(arb_radref(x));
}

ARB_INLINE void
arb_zero_pm_one(arb_t x)
{
    arf_zero(arb_midref(x));
    mag_one(arb_radref(x));
}

ARB_INLINE void
arb_unit_interval(arb_t x)
{
    arf_one(arb_midref(x));
    mag_one(arb_radref(x));
    ARF_EXP(arb_midref(x))--;
    MAG_EXP(arb_radref(x))--;
}

void arb_indeterminate(arb_t x);

int arb_is_finite(const arb_t x);

void arb_set(arb_t x, const arb_t y);

ARB_INLINE void
arb_swap(arb_t x, arb_t y)
{
    FLINT_SWAP(arb_struct, *x, *y);
}

void arb_set_round(arb_t z, const arb_t x, slong prec);

void arb_trim(arb_t y, const arb_t x);

void arb_neg(arb_t y, const arb_t x);

void arb_neg_round(arb_t x, const arb_t y, slong prec);

void arb_abs(arb_t y, const arb_t x);
void arb_nonnegative_abs(arb_t y, const arb_t x);

void arb_sgn(arb_t res, const arb_t x);
int arb_sgn_nonzero(const arb_t x);

void _arb_digits_round_inplace(char * s, flint_bitcnt_t * shift, fmpz_t error, slong n, arf_rnd_t rnd);

int arb_set_str(arb_t res, const char * inp, slong prec);

#define ARB_STR_MORE UWORD(1)
#define ARB_STR_NO_RADIUS UWORD(2)
#define ARB_STR_CONDENSE UWORD(16)

char * arb_get_str(const arb_t x, slong n, ulong flags);


ARB_INLINE void
arb_set_arf(arb_t x, const arf_t y)
{
    arf_set(arb_midref(x), y);
    mag_zero(arb_radref(x));
}

void arb_set_si(arb_t x, slong y);

void arb_set_ui(arb_t x, ulong y);

void arb_set_d(arb_t x, double y);

void arb_set_fmpz(arb_t x, const fmpz_t y);

ARB_INLINE void
arb_set_fmpz_2exp(arb_t x, const fmpz_t y, const fmpz_t exp)
{
    arf_set_fmpz_2exp(arb_midref(x), y, exp);
    mag_zero(arb_radref(x));
}

void arb_set_round_fmpz_2exp(arb_t y, const fmpz_t x, const fmpz_t exp, slong prec);

void arb_set_round_fmpz(arb_t y, const fmpz_t x, slong prec);

ARB_INLINE int
arb_is_one(const arb_t f)
{
    return arf_is_one(arb_midref(f)) && mag_is_zero(arb_radref(f));
}

void arb_one(arb_t f);

#ifdef FLINT_HAVE_FILE
void arb_fprint(FILE * file, const arb_t x);
void arb_fprintd(FILE * file, const arb_t x, slong digits);
void arb_fprintn(FILE * file, const arb_t x, slong digits, ulong flags);
#endif

void arb_print(const arb_t x);
void arb_printd(const arb_t x, slong digits);
void arb_printn(const arb_t x, slong digits, ulong flags);

void _arb_vec_printn(arb_srcptr vec, slong len, slong ndigits, ulong flags);
void _arb_vec_printd(arb_srcptr vec, slong len, slong ndigits);

void arb_mul_2exp_si(arb_t y, const arb_t x, slong e);

ARB_INLINE void
arb_mul_2exp_fmpz(arb_t y, const arb_t x, const fmpz_t e)
{
    arf_mul_2exp_fmpz(arb_midref(y), arb_midref(x), e);
    mag_mul_2exp_fmpz(arb_radref(y), arb_radref(x), e);
}

ARB_INLINE int
arb_is_int(const arb_t x)
{
    return mag_is_zero(arb_radref(x)) &&
           arf_is_int(arb_midref(x));
}

ARB_INLINE int
arb_is_int_2exp_si(const arb_t x, slong e)
{
    return mag_is_zero(arb_radref(x)) &&
           arf_is_int_2exp_si(arb_midref(x), e);
}

/* implementations are in arb/richcmp.c */
int arb_contains_zero(const arb_t x);
int arb_is_nonzero(const arb_t x);
int arb_is_positive(const arb_t x);
int arb_is_nonnegative(const arb_t x);
int arb_is_negative(const arb_t x);
int arb_is_nonpositive(const arb_t x);
int arb_contains_negative(const arb_t x);
int arb_contains_nonpositive(const arb_t x);
int arb_contains_positive(const arb_t x);
int arb_contains_nonnegative(const arb_t x);

void arb_get_mag_lower(mag_t z, const arb_t x);

void arb_get_mag_lower_nonnegative(mag_t z, const arb_t x);

ARB_INLINE void
arb_get_mag(mag_t z, const arb_t x)
{
    mag_t t;
    mag_init_set_arf(t, arb_midref(x));
    mag_add(z, t, arb_radref(x));
    mag_clear(t);
}

ARB_INLINE void
arb_get_mid_arb(arb_t z, const arb_t x)
{
    arf_set(arb_midref(z), arb_midref(x));
    mag_zero(arb_radref(z));
}

ARB_INLINE void
arb_get_rad_arb(arb_t z, const arb_t x)
{
    arf_set_mag(arb_midref(z), arb_radref(x));
    mag_zero(arb_radref(z));
}

void arb_get_abs_ubound_arf(arf_t u, const arb_t x, slong prec);
void arb_get_abs_lbound_arf(arf_t u, const arb_t x, slong prec);
void arb_get_ubound_arf(arf_t u, const arb_t x, slong prec);
void arb_get_lbound_arf(arf_t u, const arb_t x, slong prec);

void arb_nonnegative_part(arb_t res, const arb_t x);

slong arb_rel_error_bits(const arb_t x);

ARB_INLINE slong
arb_rel_accuracy_bits(const arb_t x)
{
    return -arb_rel_error_bits(x);
}

slong arb_rel_one_accuracy_bits(const arb_t x);

ARB_INLINE slong
arb_bits(const arb_t x)
{
    return arf_bits(arb_midref(x));
}

void arb_randtest_exact(arb_t x, flint_rand_t state, slong prec, slong mag_bits);

void arb_randtest_wide(arb_t x, flint_rand_t state, slong prec, slong mag_bits);

void arb_randtest_precise(arb_t x, flint_rand_t state, slong prec, slong mag_bits);

void arb_randtest_positive(arb_t x, flint_rand_t state, slong prec, slong mag_bits);

void arb_randtest(arb_t x, flint_rand_t state, slong prec, slong mag_bits);

void arb_randtest_special(arb_t x, flint_rand_t state, slong prec, slong mag_bits);

void arb_urandom(arb_t x, flint_rand_t state, slong prec);

void arb_add_error_arf(arb_t x, const arf_t err);

void arb_add_error_2exp_si(arb_t x, slong err);

void arb_add_error_2exp_fmpz(arb_t x, const fmpz_t err);

void arb_add_error(arb_t x, const arb_t error);

ARB_INLINE void
arb_add_error_mag(arb_t x, const mag_t err)
{
    mag_add(arb_radref(x), arb_radref(x), err);
}

int arb_contains_arf(const arb_t x, const arf_t y);

int arb_contains_fmpq(const arb_t x, const fmpq_t y);

int arb_contains_fmpz(const arb_t x, const fmpz_t y);

int arb_contains_si(const arb_t x, slong y);

#ifdef __MPFR_H
int arb_contains_mpfr(const arb_t x, const mpfr_t y);
#endif

int arb_overlaps(const arb_t x, const arb_t y);

int arb_contains(const arb_t x, const arb_t y);

int arb_contains_interior(const arb_t x, const arb_t y);

int arb_contains_int(const arb_t x);

void arb_get_interval_fmpz_2exp(fmpz_t a, fmpz_t b, fmpz_t exp, const arb_t x);
int arb_get_unique_fmpz(fmpz_t z, const arb_t x);

void arb_get_fmpz_mid_rad_10exp(fmpz_t mid, fmpz_t rad, fmpz_t exp, const arb_t x, slong n);

void arb_floor(arb_t z, const arb_t x, slong prec);
void arb_ceil(arb_t z, const arb_t x, slong prec);
void arb_nint(arb_t res, const arb_t x, slong prec);
void arb_trunc(arb_t res, const arb_t x, slong prec);

void arb_set_interval_arf(arb_t x, const arf_t a, const arf_t b, slong prec);
#ifdef __MPFR_H
void arb_set_interval_mpfr(arb_t x, const mpfr_t a, const mpfr_t b, slong prec);
#endif

void arb_get_interval_arf(arf_t a, arf_t b, const arb_t x, slong prec);
#ifdef __MPFR_H
void arb_get_interval_mpfr(mpfr_t a, mpfr_t b, const arb_t x);
#endif

void arb_set_interval_mag(arb_t res, const mag_t a, const mag_t b, slong prec);
void arb_set_interval_neg_pos_mag(arb_t res, const mag_t a, const mag_t b, slong prec);

void arb_union(arb_t z, const arb_t x, const arb_t y, slong prec);
int arb_intersection(arb_t z, const arb_t x, const arb_t y, slong prec);
void arb_get_rand_fmpq(fmpq_t q, flint_rand_t state, const arb_t x, slong bits);

void arb_min(arb_t z, const arb_t x, const arb_t y, slong prec);
void arb_max(arb_t z, const arb_t x, const arb_t y, slong prec);
void arb_minmax(arb_t z1, arb_t z2, const arb_t x, const arb_t y, slong prec);

int arb_can_round_arf(const arb_t x, slong prec, arf_rnd_t rnd);
#ifdef __MPFR_H
int arb_can_round_mpfr(const arb_t x, slong prec, mpfr_rnd_t rnd);
#endif

void arb_add(arb_t z, const arb_t x, const arb_t y, slong prec);
void arb_add_arf(arb_t z, const arb_t x, const arf_t y, slong prec);
void arb_add_ui(arb_t z, const arb_t x, ulong y, slong prec);
void arb_add_si(arb_t z, const arb_t x, slong y, slong prec);
void arb_add_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec);
void arb_add_fmpz_2exp(arb_t z, const arb_t x, const fmpz_t man, const fmpz_t exp, slong prec);

void arb_sub(arb_t z, const arb_t x, const arb_t y, slong prec);
void arb_sub_arf(arb_t z, const arb_t x, const arf_t y, slong prec);
void arb_sub_ui(arb_t z, const arb_t x, ulong y, slong prec);
void arb_sub_si(arb_t z, const arb_t x, slong y, slong prec);
void arb_sub_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec);

void arb_mul(arb_t z, const arb_t x, const arb_t y, slong prec);
void arb_mul_arf(arb_t z, const arb_t x, const arf_t y, slong prec);
void arb_mul_si(arb_t z, const arb_t x, slong y, slong prec);
void arb_mul_ui(arb_t z, const arb_t x, ulong y, slong prec);
void arb_mul_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec);

void arb_addmul(arb_t z, const arb_t x, const arb_t y, slong prec);
void arb_addmul_arf(arb_t z, const arb_t x, const arf_t y, slong prec);
void arb_addmul_si(arb_t z, const arb_t x, slong y, slong prec);
void arb_addmul_ui(arb_t z, const arb_t x, ulong y, slong prec);
void arb_addmul_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec);

void arb_submul(arb_t z, const arb_t x, const arb_t y, slong prec);
void arb_submul_arf(arb_t z, const arb_t x, const arf_t y, slong prec);
void arb_submul_si(arb_t z, const arb_t x, slong y, slong prec);
void arb_submul_ui(arb_t z, const arb_t x, ulong y, slong prec);
void arb_submul_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec);

void arb_fma(arb_t res, const arb_t x, const arb_t y, const arb_t z, slong prec);
void arb_fma_arf(arb_t res, const arb_t x, const arf_t y, const arb_t z, slong prec);
void arb_fma_si(arb_t res, const arb_t x, slong y, const arb_t z, slong prec);
void arb_fma_ui(arb_t res, const arb_t x, ulong y, const arb_t z, slong prec);
void arb_fma_fmpz(arb_t res, const arb_t x, const fmpz_t y, const arb_t z, slong prec);

void arb_dot_simple(arb_t res, const arb_t initial, int subtract,
    arb_srcptr x, slong xstep, arb_srcptr y, slong ystep, slong len, slong prec);
void arb_dot_precise(arb_t res, const arb_t initial, int subtract,
    arb_srcptr x, slong xstep, arb_srcptr y, slong ystep, slong len, slong prec);
void arb_dot(arb_t res, const arb_t initial, int subtract,
    arb_srcptr x, slong xstep, arb_srcptr y, slong ystep, slong len, slong prec);

void arb_approx_dot(arb_t res, const arb_t initial, int subtract,
    arb_srcptr x, slong xstep, arb_srcptr y, slong ystep, slong len, slong prec);

void arb_dot_ui(arb_t res, const arb_t initial, int subtract,
    arb_srcptr x, slong xstep, const ulong * y, slong ystep, slong len, slong prec);
void arb_dot_si(arb_t res, const arb_t initial, int subtract,
    arb_srcptr x, slong xstep, const slong * y, slong ystep, slong len, slong prec);
void arb_dot_uiui(arb_t res, const arb_t initial, int subtract,
    arb_srcptr x, slong xstep, const ulong * y, slong ystep, slong len, slong prec);
void arb_dot_siui(arb_t res, const arb_t initial, int subtract,
    arb_srcptr x, slong xstep, const ulong * y, slong ystep, slong len, slong prec);
void arb_dot_fmpz(arb_t res, const arb_t initial, int subtract,
    arb_srcptr x, slong xstep, const fmpz * y, slong ystep, slong len, slong prec);

void arb_div(arb_t z, const arb_t x, const arb_t y, slong prec);
void arb_div_arf(arb_t z, const arb_t x, const arf_t y, slong prec);
void arb_div_si(arb_t z, const arb_t x, slong y, slong prec);
void arb_div_ui(arb_t z, const arb_t x, ulong y, slong prec);
void arb_div_fmpz(arb_t z, const arb_t x, const fmpz_t y, slong prec);
void arb_fmpz_div_fmpz(arb_t z, const fmpz_t x, const fmpz_t y, slong prec);
void arb_ui_div(arb_t z, ulong x, const arb_t y, slong prec);

void arb_div_arf_newton(arb_t res, const arb_t x, const arf_t y, slong prec);
void arb_div_newton(arb_t res, const arb_t x, const arb_t y, slong prec);
void arb_fmpz_divapprox(fmpz_t res, const fmpz_t x, const fmpz_t y);

ARB_INLINE void
arb_inv(arb_t y, const arb_t x, slong prec)
{
    arb_ui_div(y, 1, x, prec);
}

ARB_INLINE void
arb_set_fmpq(arb_t y, const fmpq_t x, slong prec)
{
    arb_fmpz_div_fmpz(y, fmpq_numref(x), fmpq_denref(x), prec);
}

void arb_sqrt(arb_t z, const arb_t x, slong prec);
void arb_sqrt_arf(arb_t z, const arf_t x, slong prec);
void arb_sqrt_fmpz(arb_t z, const fmpz_t x, slong prec);
void arb_sqrt_ui(arb_t z, ulong x, slong prec);

void arb_sqrtpos(arb_t z, const arb_t x, slong prec);
void arb_hypot(arb_t z, const arb_t x, const arb_t y, slong prec);

void arb_rsqrt(arb_t z, const arb_t x, slong prec);
void arb_rsqrt_ui(arb_t z, ulong x, slong prec);
void arb_sqrt1pm1(arb_t r, const arb_t z, slong prec);

void arb_rsqrt_arf_newton(arb_t res, const arf_t x, slong prec);
void arb_sqrt_arf_newton(arb_t res, const arf_t x, slong prec);
void arb_sqrt_newton(arb_t z, const arb_t x, slong prec);

void arb_pow_fmpz_binexp(arb_t y, const arb_t b, const fmpz_t e, slong prec);
void arb_pow_fmpz(arb_t y, const arb_t b, const fmpz_t e, slong prec);
void arb_pow_ui(arb_t y, const arb_t b, ulong e, slong prec);
void arb_ui_pow_ui(arb_t y, ulong b, ulong e, slong prec);
void arb_si_pow_ui(arb_t y, slong b, ulong e, slong prec);
void arb_pow_fmpq(arb_t y, const arb_t x, const fmpq_t a, slong prec);

void arb_div_2expm1_ui(arb_t z, const arb_t x, ulong n, slong prec);
void arb_pow(arb_t z, const arb_t x, const arb_t y, slong prec);
void arb_root_ui(arb_t z, const arb_t x, ulong k, slong prec);
void arb_root(arb_t z, const arb_t x, ulong k, slong prec); /* back compat */
void arb_log(arb_t z, const arb_t x, slong prec);
void arb_log_arf(arb_t z, const arf_t x, slong prec);
void arb_log_ui(arb_t z, ulong x, slong prec);
void arb_log_fmpz(arb_t z, const fmpz_t x, slong prec);
void arb_log1p(arb_t r, const arb_t z, slong prec);
void arb_log_base_ui(arb_t res, const arb_t x, ulong b, slong prec);
void arb_log_hypot(arb_t res, const arb_t x, const arb_t y, slong prec);
void arb_exp(arb_t z, const arb_t x, slong prec);
void arb_expm1(arb_t z, const arb_t x, slong prec);
void arb_exp_invexp(arb_t z, arb_t w, const arb_t x, slong prec);
void arb_sin(arb_t s, const arb_t x, slong prec);
void arb_cos(arb_t c, const arb_t x, slong prec);
void arb_sin_cos(arb_t s, arb_t c, const arb_t x, slong prec);
void arb_sin_pi(arb_t s, const arb_t x, slong prec);
void arb_cos_pi(arb_t c, const arb_t x, slong prec);
void arb_sin_cos_pi(arb_t s, arb_t c, const arb_t x, slong prec);
void arb_tan(arb_t y, const arb_t x, slong prec);
void arb_cot(arb_t y, const arb_t x, slong prec);
void arb_tan_pi(arb_t y, const arb_t x, slong prec);
void arb_cot_pi(arb_t y, const arb_t x, slong prec);
void _arb_sin_pi_fmpq_algebraic(arb_t s, ulong p, ulong q, slong prec);
void _arb_cos_pi_fmpq_algebraic(arb_t c, ulong p, ulong q, slong prec);
void _arb_sin_cos_pi_fmpq_algebraic(arb_t s, arb_t c, ulong p, ulong q, slong prec);
void arb_sin_cos_pi_fmpq(arb_t s, arb_t c, const fmpq_t x, slong prec);
void arb_sin_pi_fmpq(arb_t s, const fmpq_t x, slong prec);
void arb_cos_pi_fmpq(arb_t c, const fmpq_t x, slong prec);
void arb_sinc(arb_t z, const arb_t x, slong prec);
void arb_sinc_pi(arb_t z, const arb_t x, slong prec);
void arb_sinh(arb_t z, const arb_t x, slong prec);
void arb_cosh(arb_t z, const arb_t x, slong prec);
void arb_sinh_cosh(arb_t s, arb_t c, const arb_t x, slong prec);
void arb_tanh(arb_t y, const arb_t x, slong prec);
void arb_coth(arb_t y, const arb_t x, slong prec);
void arb_atan_arf(arb_t z, const arf_t x, slong prec);
void arb_atan(arb_t z, const arb_t x, slong prec);
void arb_atan2(arb_t z, const arb_t b, const arb_t a, slong prec);
void arb_asin(arb_t z, const arb_t x, slong prec);
void arb_acos(arb_t z, const arb_t x, slong prec);
void arb_atanh(arb_t z, const arb_t x, slong prec);
void arb_asinh(arb_t z, const arb_t x, slong prec);
void arb_acosh(arb_t z, const arb_t x, slong prec);

void arb_sec(arb_t res, const arb_t x, slong prec);
void arb_csc(arb_t res, const arb_t x, slong prec);
void arb_csc_pi(arb_t res, const arb_t x, slong prec);
void arb_sech(arb_t res, const arb_t x, slong prec);
void arb_csch(arb_t res, const arb_t x, slong prec);

void arb_fac_ui(arb_t z, ulong n, slong prec);
void arb_doublefac_ui(arb_t z, ulong n, slong prec);
void arb_bin_ui(arb_t z, const arb_t n, ulong k, slong prec);
void arb_bin_uiui(arb_t z, ulong n, ulong k, slong prec);
void arb_fib_fmpz(arb_t z, const fmpz_t n, slong prec);
void arb_fib_ui(arb_t z, ulong n, slong prec);
void arb_const_pi(arb_t z, slong prec);
void arb_const_sqrt_pi(arb_t z, slong prec);
void arb_const_log_sqrt2pi(arb_t z, slong prec);
void arb_const_log2(arb_t z, slong prec);
void arb_const_log10(arb_t z, slong prec);
void arb_const_euler(arb_t z, slong prec);
void arb_const_catalan(arb_t z, slong prec);
void arb_const_e(arb_t z, slong prec);
void arb_const_khinchin(arb_t z, slong prec);
void arb_const_glaisher(arb_t z, slong prec);
void arb_const_reciprocal_fibonacci(arb_t res, slong prec);
void arb_agm(arb_t z, const arb_t x, const arb_t y, slong prec);
void arb_lgamma(arb_t z, const arb_t x, slong prec);
void arb_rgamma(arb_t z, const arb_t x, slong prec);
void arb_gamma(arb_t z, const arb_t x, slong prec);
void arb_gamma_fmpq(arb_t z, const fmpq_t x, slong prec);
void arb_gamma_fmpz(arb_t z, const fmpz_t x, slong prec);
void arb_digamma(arb_t y, const arb_t x, slong prec);
void arb_zeta(arb_t z, const arb_t s, slong prec);
void arb_hurwitz_zeta(arb_t z, const arb_t s, const arb_t a, slong prec);

void arb_rising_ui(arb_t z, const arb_t x, ulong n, slong prec);
void arb_rising_fmpq_ui(arb_t y, const fmpq_t x, ulong n, slong prec);
void arb_rising(arb_t z, const arb_t x, const arb_t n, slong prec);
void arb_rising2_ui(arb_t u, arb_t v, const arb_t x, ulong n, slong prec);

void arb_log_ui_from_prev(arb_t s, ulong k, arb_t log_prev, ulong prev, slong prec);

void arb_const_apery(arb_t s, slong prec);

void arb_zeta_ui_asymp(arb_t x, ulong s, slong prec);
void arb_zeta_ui_borwein_bsplit(arb_t x, ulong s, slong prec);
void arb_zeta_ui_euler_product(arb_t z, ulong s, slong prec);
void arb_zeta_ui_bernoulli(arb_t x, ulong n, slong prec);
void arb_zeta_ui_vec_borwein(arb_ptr z, ulong start, slong num, ulong step, slong prec);
void arb_zeta_ui(arb_t x, ulong n, slong prec);
void arb_zeta_ui_vec_even(arb_ptr x, ulong start, slong num, slong prec);
void arb_zeta_ui_vec_odd(arb_ptr x, ulong start, slong num, slong prec);
void arb_zeta_ui_vec(arb_ptr x, ulong start, slong num, slong prec);
void arb_bernoulli_ui(arb_t b, ulong n, slong prec);
void arb_bernoulli_ui_zeta(arb_t b, ulong n, slong prec);
void arb_bernoulli_fmpz(arb_t b, const fmpz_t n, slong prec);

void arb_bernoulli_poly_ui(arb_t res, ulong n, const arb_t x, slong prec);

void arb_polylog(arb_t w, const arb_t s, const arb_t z, slong prec);
void arb_polylog_si(arb_t w, slong s, const arb_t z, slong prec);

void arb_chebyshev_t_ui(arb_t a, ulong n, const arb_t x, slong prec);
void arb_chebyshev_t2_ui(arb_t a, arb_t b, ulong n, const arb_t x, slong prec);
void arb_chebyshev_u_ui(arb_t a, ulong n, const arb_t x, slong prec);
void arb_chebyshev_u2_ui(arb_t a, arb_t b, ulong n, const arb_t x, slong prec);

void arb_power_sum_vec(arb_ptr res, const arb_t a, const arb_t b, slong len, slong prec);
void arb_bell_sum_taylor(arb_t res, const fmpz_t n, const fmpz_t a, const fmpz_t b, const fmpz_t mmag, slong prec);
void arb_bell_sum_bsplit(arb_t res, const fmpz_t n, const fmpz_t a, const fmpz_t b, const fmpz_t mmag, slong prec);
void arb_bell_fmpz(arb_t res, const fmpz_t n, slong prec);
void arb_bell_ui(arb_t res, ulong n, slong prec);

void arb_euler_number_fmpz(arb_t res, const fmpz_t n, slong prec);
void arb_euler_number_ui(arb_t res, ulong n, slong prec);
void arb_fmpz_euler_number_ui_multi_mod(fmpz_t res, ulong n, double alpha);
void arb_fmpz_euler_number_ui(fmpz_t res, ulong n);

void arb_partitions_fmpz(arb_t res, const fmpz_t n, slong prec);
void arb_partitions_ui(arb_t res, ulong n, slong prec);

void arb_primorial_nth_ui(arb_t res, ulong n, slong prec);
void arb_primorial_ui(arb_t res, ulong n, slong prec);

void arb_lambertw(arb_t res, const arb_t x, int flags, slong prec);

ARB_INLINE void
arb_sqr(arb_t res, const arb_t val, slong prec)
{
    arb_mul(res, val, val, prec);
}

#define ARB_DEF_CACHED_CONSTANT(name, comp_func) \
    FLINT_TLS_PREFIX slong name ## _cached_prec = 0; \
    FLINT_TLS_PREFIX arb_t name ## _cached_value; \
    void name ## _cleanup(void) \
    { \
        arb_clear(name ## _cached_value); \
        name ## _cached_prec = 0; \
    } \
    void name(arb_t x, slong prec) \
    { \
        if (name ## _cached_prec < prec) \
        { \
            if (name ## _cached_prec == 0) \
            { \
                arb_init(name ## _cached_value); \
                flint_register_cleanup_function(name ## _cleanup); \
            } \
            comp_func(name ## _cached_value, prec + 32); \
            name ## _cached_prec = prec; \
        } \
        arb_set_round(x, name ## _cached_value, prec); \
    }

/* vector functions */

ARB_INLINE arb_ptr
_arb_vec_entry_ptr(arb_ptr vec, slong i)
{
    return vec + i;
}

ARB_INLINE void
_arb_vec_zero(arb_ptr A, slong n)
{
    slong i;
    for (i = 0; i < n; i++)
        arb_zero(A + i);
}

ARB_INLINE int
_arb_vec_is_zero(arb_srcptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        if (!arb_is_zero(vec + i))
            return 0;
    return 1;
}

ARB_INLINE int
_arb_vec_is_finite(arb_srcptr x, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
        if (!arb_is_finite(x + i))
            return 0;

    return 1;
}

ARB_INLINE int
_arb_vec_equal(arb_srcptr vec1, arb_srcptr vec2, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (!arb_equal(vec1 + i, vec2 + i))
            return 0;
    }
    return 1;
}

ARB_INLINE int
_arb_vec_overlaps(arb_srcptr vec1, arb_srcptr vec2, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (!arb_overlaps(vec1 + i, vec2 + i))
            return 0;
    }

    return 1;
}

ARB_INLINE int
_arb_vec_contains(arb_srcptr vec1, arb_srcptr vec2, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (!arb_contains(vec1 + i, vec2 + i))
            return 0;
    }

    return 1;
}

ARB_INLINE void
_arb_vec_set(arb_ptr res, arb_srcptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        arb_set(res + i, vec + i);
}

ARB_INLINE void
_arb_vec_set_round(arb_ptr res, arb_srcptr vec, slong len, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        arb_set_round(res + i, vec + i, prec);
}

ARB_INLINE void
_arb_vec_swap(arb_ptr res, arb_ptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        arb_swap(res + i, vec + i);
}

ARB_INLINE void
_arb_vec_neg(arb_ptr B, arb_srcptr A, slong n)
{
    slong i;
    for (i = 0; i < n; i++)
        arb_neg(B + i, A + i);
}

ARB_INLINE void
_arb_vec_sub(arb_ptr C, arb_srcptr A,
    arb_srcptr B, slong n, slong prec)
{
    slong i;
    for (i = 0; i < n; i++)
        arb_sub(C + i, A + i, B + i, prec);
}

ARB_INLINE void
_arb_vec_add(arb_ptr C, arb_srcptr A,
    arb_srcptr B, slong n, slong prec)
{
    slong i;
    for (i = 0; i < n; i++)
        arb_add(C + i, A + i, B + i, prec);
}

ARB_INLINE void
_arb_vec_scalar_mul(arb_ptr res, arb_srcptr vec,
    slong len, const arb_t c, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        arb_mul(res + i, vec + i, c, prec);
}

ARB_INLINE void
_arb_vec_scalar_div(arb_ptr res, arb_srcptr vec,
    slong len, const arb_t c, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        arb_div(res + i, vec + i, c, prec);
}

ARB_INLINE void
_arb_vec_scalar_mul_fmpz(arb_ptr res, arb_srcptr vec,
    slong len, const fmpz_t c, slong prec)
{
    slong i;
    arf_t t;
    arf_init(t);
    arf_set_fmpz(t, c);
    for (i = 0; i < len; i++)
        arb_mul_arf(res + i, vec + i, t, prec);
    arf_clear(t);
}

ARB_INLINE void
_arb_vec_scalar_mul_2exp_si(arb_ptr res, arb_srcptr src, slong len, slong c)
{
    slong i;
    for (i = 0; i < len; i++)
        arb_mul_2exp_si(res + i, src + i, c);
}

ARB_INLINE void
_arb_vec_scalar_addmul(arb_ptr res, arb_srcptr vec,
    slong len, const arb_t c, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        arb_addmul(res + i, vec + i, c, prec);
}

void _arb_vec_get_mag(mag_t bound, arb_srcptr vec, slong len);

ARB_INLINE slong
_arb_vec_bits(arb_srcptr x, slong len)
{
    slong i, b, c;

    b = 0;
    for (i = 0; i < len; i++)
    {
        c = arb_bits(x + i);
        b = FLINT_MAX(b, c);
    }

    return b;
}

void _arb_vec_set_powers(arb_ptr xs, const arb_t x, slong len, slong prec);

ARB_INLINE void
_arb_vec_add_error_arf_vec(arb_ptr res, arf_srcptr err, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        arb_add_error_arf(res + i, err + i);
}

ARB_INLINE void
_arb_vec_add_error_mag_vec(arb_ptr res, mag_srcptr err, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        mag_add(arb_radref(res + i), arb_radref(res + i), err + i);
}

ARB_INLINE void
_arb_vec_indeterminate(arb_ptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        arb_indeterminate(vec + i);
}

ARB_INLINE void
_arb_vec_trim(arb_ptr res, arb_srcptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        arb_trim(res + i, vec + i);
}

ARB_INLINE int
_arb_vec_get_unique_fmpz_vec(fmpz * res,  arb_srcptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        if (!arb_get_unique_fmpz(res + i, vec + i))
            return 0;
    return 1;
}

/* arctangent implementation */

#define ARB_ATAN_TAB1_BITS 8
#define ARB_ATAN_TAB1_PREC 512
#define ARB_ATAN_TAB1_LIMBS (ARB_ATAN_TAB1_PREC / FLINT_BITS)

#define ARB_ATAN_TAB21_BITS 5
#define ARB_ATAN_TAB22_BITS 5
#define ARB_ATAN_TAB2_PREC 4608
#define ARB_ATAN_TAB2_LIMBS (ARB_ATAN_TAB2_PREC / FLINT_BITS)

FLINT_DLL extern const mp_limb_t arb_atan_tab1[1 << ARB_ATAN_TAB1_BITS][ARB_ATAN_TAB1_LIMBS];
FLINT_DLL extern const mp_limb_t arb_atan_tab21[1 << ARB_ATAN_TAB21_BITS][ARB_ATAN_TAB2_LIMBS];
FLINT_DLL extern const mp_limb_t arb_atan_tab22[1 << ARB_ATAN_TAB22_BITS][ARB_ATAN_TAB2_LIMBS];
FLINT_DLL extern const mp_limb_t arb_atan_pi2_minus_one[ARB_ATAN_TAB2_LIMBS];

void
_arb_atan_taylor_naive(mp_ptr y, mp_limb_t * error,
    mp_srcptr x, mp_size_t xn, ulong N, int alternating);

void _arb_atan_taylor_rs(mp_ptr y, mp_limb_t * error,
    mp_srcptr x, mp_size_t xn, ulong N, int alternating);

#define ARB_ATAN_NEWTON_PREC 3400

void arb_atan_newton(arb_t res, const arb_t x, slong prec);
void arb_atan_arf_newton(arb_t res, const arf_t x, slong prec);

/* logarithm implementation */

#define ARB_LOG_TAB11_BITS 7
#define ARB_LOG_TAB12_BITS 7
#define ARB_LOG_TAB1_PREC 512
#define ARB_LOG_TAB1_LIMBS (ARB_LOG_TAB1_PREC / FLINT_BITS)

#define ARB_LOG_TAB21_BITS 5
#define ARB_LOG_TAB22_BITS 5
#define ARB_LOG_TAB2_PREC 4608
#define ARB_LOG_TAB2_LIMBS (ARB_LOG_TAB2_PREC / FLINT_BITS)

FLINT_DLL extern const mp_limb_t arb_log_tab11[1 << ARB_LOG_TAB11_BITS][ARB_LOG_TAB1_LIMBS];
FLINT_DLL extern const mp_limb_t arb_log_tab12[1 << ARB_LOG_TAB12_BITS][ARB_LOG_TAB1_LIMBS];
FLINT_DLL extern const mp_limb_t arb_log_tab21[1 << ARB_LOG_TAB21_BITS][ARB_LOG_TAB2_LIMBS];
FLINT_DLL extern const mp_limb_t arb_log_tab22[1 << ARB_LOG_TAB22_BITS][ARB_LOG_TAB2_LIMBS];
FLINT_DLL extern const mp_srcptr arb_log_log2_tab;

void arb_log_newton(arb_t res, const arb_t x, slong prec);
void arb_log_arf_newton(arb_t res, const arf_t x, slong prec);

/* logarithms of primes */

#define ARB_LOG_PRIME_CACHE_NUM 13

FLINT_DLL extern const mp_limb_t arb_log_p_tab[ARB_LOG_PRIME_CACHE_NUM][ARB_LOG_TAB2_LIMBS];
void arb_log_primes_vec_bsplit(arb_ptr res, slong n, slong prec);

void _arb_log_p_ensure_cached(slong prec);
arb_srcptr _arb_log_p_cache_vec(void);

/* exponential implementation */

/* only goes up to log(2) * 256 */
#define ARB_EXP_TAB1_NUM 178
#define ARB_EXP_TAB1_BITS 8
#define ARB_EXP_TAB1_PREC 512
#define ARB_EXP_TAB1_LIMBS (ARB_EXP_TAB1_PREC / FLINT_BITS)

/* only goes up to log(2) * 32 */
#define ARB_EXP_TAB21_NUM 23
#define ARB_EXP_TAB21_BITS 5
#define ARB_EXP_TAB22_NUM (1 << ARB_EXP_TAB22_BITS)
#define ARB_EXP_TAB22_BITS 5
#define ARB_EXP_TAB2_PREC 4608
#define ARB_EXP_TAB2_LIMBS (ARB_EXP_TAB2_PREC / FLINT_BITS)

FLINT_DLL extern const mp_limb_t arb_exp_tab1[ARB_EXP_TAB1_NUM][ARB_EXP_TAB1_LIMBS];
FLINT_DLL extern const mp_limb_t arb_exp_tab21[ARB_EXP_TAB21_NUM][ARB_EXP_TAB2_LIMBS];
FLINT_DLL extern const mp_limb_t arb_exp_tab22[ARB_EXP_TAB22_NUM][ARB_EXP_TAB2_LIMBS];

void _arb_exp_taylor_naive(mp_ptr y, mp_limb_t * error,
    mp_srcptr x, mp_size_t xn, ulong N);

void _arb_exp_taylor_rs(mp_ptr y, mp_limb_t * error,
    mp_srcptr x, mp_size_t xn, ulong N);

void arb_exp_arf_bb(arb_t z, const arf_t x, slong prec, int minus_one);
void arb_exp_arf_rs_generic(arb_t res, const arf_t x, slong prec, int minus_one);

int _arb_get_mpn_fixed_mod_log2(mp_ptr w, fmpz_t q, mp_limb_t * error,
                                                const arf_t x, mp_size_t wn);

slong _arb_exp_taylor_bound(slong mag, slong prec);

void _arb_exp_sum_bs_powtab(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp,
    const fmpz_t x, flint_bitcnt_t r, slong N);

void _arb_exp_sum_bs_simple(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp,
    const fmpz_t x, flint_bitcnt_t r, slong N);

#define ARB_LOG_REDUCTION_DEFAULT_MAX_PREC 4000000
#define ARB_EXP_LOG_REDUCTION_PREC 2240
#define ARB_LOG_NEWTON_PREC 2800

void arb_exp_arf_log_reduction(arb_t res, const arf_t x, slong prec, int minus_one);
void arb_exp_arf_generic(arb_t z, const arf_t x, slong prec, int minus_one);
void arb_exp_arf(arb_t z, const arf_t x, slong prec, int minus_one, slong maglim);

/* sin/cos implementation */

/* only goes up to (pi/4) * 256 */
#define ARB_SIN_COS_TAB1_NUM 203
#define ARB_SIN_COS_TAB1_BITS 8
#define ARB_SIN_COS_TAB1_PREC 512
#define ARB_SIN_COS_TAB1_LIMBS (ARB_SIN_COS_TAB1_PREC / FLINT_BITS)

/* only goes up to (pi/4) * 32 */
#define ARB_SIN_COS_TAB21_NUM 26
#define ARB_SIN_COS_TAB21_BITS 5
#define ARB_SIN_COS_TAB22_NUM (1 << ARB_SIN_COS_TAB22_BITS)
#define ARB_SIN_COS_TAB22_BITS 5
#define ARB_SIN_COS_TAB2_PREC 4608
#define ARB_SIN_COS_TAB2_LIMBS (ARB_SIN_COS_TAB2_PREC / FLINT_BITS)

FLINT_DLL extern const mp_limb_t arb_sin_cos_tab1[2 * ARB_SIN_COS_TAB1_NUM][ARB_SIN_COS_TAB1_LIMBS];
FLINT_DLL extern const mp_limb_t arb_sin_cos_tab21[2 * ARB_SIN_COS_TAB21_NUM][ARB_SIN_COS_TAB2_LIMBS];
FLINT_DLL extern const mp_limb_t arb_sin_cos_tab22[2 * ARB_SIN_COS_TAB22_NUM][ARB_SIN_COS_TAB2_LIMBS];

#define ARB_PI4_TAB_LIMBS (4608 / FLINT_BITS)
FLINT_DLL extern const mp_limb_t arb_pi4_tab[ARB_PI4_TAB_LIMBS];

void _arb_sin_cos_taylor_naive(mp_ptr ysin, mp_ptr ycos, mp_limb_t * error,
    mp_srcptr x, mp_size_t xn, ulong N);

void _arb_sin_cos_taylor_rs(mp_ptr ysin, mp_ptr ycos,
    mp_limb_t * error, mp_srcptr x, mp_size_t xn, ulong N,
    int sinonly, int alternating);

int _arb_get_mpn_fixed_mod_pi4(mp_ptr w, fmpz_t q, int * octant,
    mp_limb_t * error, const arf_t x, mp_size_t wn);

void arb_sin_cos_arf_bb(arb_t zsin, arb_t zcos, const arf_t x, slong prec);
void arb_sin_cos_arf_rs_generic(arb_t res_sin, arb_t res_cos, const arf_t x, slong prec);
void arb_sin_cos_arf_generic(arb_t res_sin, arb_t res_cos, const arf_t x, slong prec);

void _arb_sin_cos_wide(arb_t s, arb_t c, const arf_t x, const mag_t r, slong prec);
void arb_sin_cos_wide(arb_t s, arb_t c, const arb_t x, slong prec);

void _arb_sin_cos_generic(arb_t s, arb_t c, const arf_t x, const mag_t xrad, slong prec);
void arb_sin_cos_generic(arb_t s, arb_t c, const arb_t x, slong prec);

void arb_atan_gauss_primes_vec_bsplit(arb_ptr res, slong n, slong prec);

#define ARB_ATAN_GAUSS_PRIME_CACHE_NUM 13

#define ARB_SIN_COS_ATAN_REDUCTION_DEFAULT_MAX_PREC 4000000
#define ARB_SIN_COS_ATAN_REDUCTION_PREC 2600

FLINT_DLL extern const mp_limb_t arb_atan_gauss_tab[ARB_ATAN_GAUSS_PRIME_CACHE_NUM][ARB_ATAN_TAB2_LIMBS];

void _arb_atan_gauss_p_ensure_cached(slong prec);
arb_srcptr _arb_atan_gauss_p_cache_vec(void);

void arb_sin_cos_arf_atan_reduction(arb_t res1, arb_t res2, const arf_t x, slong prec);


ARB_INLINE flint_bitcnt_t
_arb_mpn_leading_zeros(mp_srcptr d, mp_size_t n)
{
    mp_limb_t t;
    mp_size_t zero_limbs;
    flint_bitcnt_t bits;

    zero_limbs = 0;

    while (1)
    {
        t = d[n - zero_limbs - 1];

        if (t != 0)
        {
            bits = flint_clz(t);
            return bits + FLINT_BITS * zero_limbs;
        }

        zero_limbs++;

        if (zero_limbs == n)
            return FLINT_BITS * zero_limbs;
    }
}

void _arb_atan_sum_bs_simple(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp,
    const fmpz_t x, flint_bitcnt_t r, slong N);

void _arb_atan_sum_bs_powtab(fmpz_t T, fmpz_t Q, flint_bitcnt_t * Qexp,
    const fmpz_t x, flint_bitcnt_t r, slong N);

void arb_atan_arf_bb(arb_t z, const arf_t x, slong prec);

void arb_atan_frac_bsplit(arb_t s, const fmpz_t p, const fmpz_t q, int hyperbolic, slong prec);

ARB_INLINE slong
arb_allocated_bytes(const arb_t x)
{
    return arf_allocated_bytes(arb_midref(x)) + mag_allocated_bytes(arb_radref(x));
}

ARB_INLINE slong
_arb_vec_allocated_bytes(arb_srcptr vec, slong len)
{
    slong i, size;

    size = len * sizeof(arb_struct);

    for (i = 0; i < len; i++)
        size += arb_allocated_bytes(vec + i);

    return size;
}

ARB_INLINE double
_arb_vec_estimate_allocated_bytes(slong len, slong prec)
{
    double size;

    size = len * (double) sizeof(arb_struct);

    if (prec > ARF_NOPTR_LIMBS * FLINT_BITS)
        size += len * (double) ((prec + FLINT_BITS - 1) / FLINT_BITS) * sizeof(mp_limb_t);

    return size;
}

int arb_load_str(arb_t res, const char * data);
char * arb_dump_str(const arb_t x);

#ifdef FLINT_HAVE_FILE
int arb_load_file(arb_t res, FILE *stream);
int arb_dump_file(FILE* stream, const arb_t x);
#endif

#ifdef __cplusplus
}
#endif

#endif
