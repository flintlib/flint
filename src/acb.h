/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_H
#define ACB_H

#ifdef ACB_INLINES_C
#define ACB_INLINE
#else
#define ACB_INLINE static inline
#endif

#include "arb.h"
#include "acb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define acb_realref(x) (&(x)->real)
#define acb_imagref(x) (&(x)->imag)

ACB_INLINE void
acb_init(acb_t x)
{
    arb_init(acb_realref(x));
    arb_init(acb_imagref(x));
}

void acb_clear(acb_t x);

acb_ptr _acb_vec_init(slong n);
void _acb_vec_clear(acb_ptr v, slong n);

ACB_INLINE arb_ptr acb_real_ptr(acb_t z) { return acb_realref(z); }
ACB_INLINE arb_ptr acb_imag_ptr(acb_t z) { return acb_imagref(z); }

ACB_INLINE void
acb_get_real(arb_t re, const acb_t z)
{
    arb_set(re, acb_realref(z));
}

ACB_INLINE void
acb_get_imag(arb_t im, const acb_t z)
{
    arb_set(im, acb_imagref(z));
}

ACB_INLINE void
acb_get_mid(acb_t res, const acb_t x)
{
    arb_get_mid_arb(acb_realref(res), acb_realref(x));
    arb_get_mid_arb(acb_imagref(res), acb_imagref(x));
}

ACB_INLINE int
acb_is_zero(const acb_t z)
{
    return arb_is_zero(acb_realref(z)) && arb_is_zero(acb_imagref(z));
}

ACB_INLINE int
acb_is_one(const acb_t z)
{
    return arb_is_one(acb_realref(z)) && arb_is_zero(acb_imagref(z));
}

ACB_INLINE int
acb_is_exact(const acb_t z)
{
    return arb_is_exact(acb_realref(z)) && arb_is_exact(acb_imagref(z));
}

ACB_INLINE int
acb_is_int(const acb_t z)
{
    return arb_is_zero(acb_imagref(z)) && arb_is_int(acb_realref(z));
}

ACB_INLINE int
acb_is_int_2exp_si(const acb_t z, slong e)
{
    return arb_is_zero(acb_imagref(z)) && arb_is_int_2exp_si(acb_realref(z), e);
}

ACB_INLINE void
acb_zero(acb_t z)
{
    arb_zero(acb_realref(z));
    arb_zero(acb_imagref(z));
}

ACB_INLINE void
acb_one(acb_t z)
{
    arb_one(acb_realref(z));
    arb_zero(acb_imagref(z));
}

ACB_INLINE void
acb_onei(acb_t z)
{
    arb_zero(acb_realref(z));
    arb_one(acb_imagref(z));
}

ACB_INLINE void
acb_set(acb_t z, const acb_t x)
{
    arb_set(acb_realref(z), acb_realref(x));
    arb_set(acb_imagref(z), acb_imagref(x));
}

ACB_INLINE void
acb_set_round(acb_t z, const acb_t x, slong prec)
{
    arb_set_round(acb_realref(z), acb_realref(x), prec);
    arb_set_round(acb_imagref(z), acb_imagref(x), prec);
}

ACB_INLINE void
acb_neg_round(acb_t z, const acb_t x, slong prec)
{
    arb_neg_round(acb_realref(z), acb_realref(x), prec);
    arb_neg_round(acb_imagref(z), acb_imagref(x), prec);
}

ACB_INLINE void
acb_swap(acb_t z, acb_t x)
{
    arb_swap(acb_realref(z), acb_realref(x));
    arb_swap(acb_imagref(z), acb_imagref(x));
}

ACB_INLINE int
acb_equal(const acb_t x, const acb_t y)
{
    return arb_equal(acb_realref(x), acb_realref(y)) &&
           arb_equal(acb_imagref(x), acb_imagref(y));
}

ACB_INLINE int
acb_equal_si(const acb_t x, slong y)
{
    return arb_equal_si(acb_realref(x), y) && arb_is_zero(acb_imagref(x));
}

ACB_INLINE int
acb_eq(const acb_t x, const acb_t y)
{
    return arb_eq(acb_realref(x), acb_realref(y)) &&
           arb_eq(acb_imagref(x), acb_imagref(y));
}

ACB_INLINE int
acb_ne(const acb_t x, const acb_t y)
{
    return arb_ne(acb_realref(x), acb_realref(y)) ||
           arb_ne(acb_imagref(x), acb_imagref(y));
}

ACB_INLINE int
acb_overlaps(const acb_t x, const acb_t y)
{
    return arb_overlaps(acb_realref(x), acb_realref(y)) &&
            arb_overlaps(acb_imagref(x), acb_imagref(y));
}

ACB_INLINE int
acb_contains_zero(const acb_t x)
{
    return arb_contains_zero(acb_realref(x)) &&
            arb_contains_zero(acb_imagref(x));
}

ACB_INLINE int
acb_contains_fmpq(const acb_t x, const fmpq_t y)
{
    return arb_contains_fmpq(acb_realref(x), y) &&
            arb_contains_zero(acb_imagref(x));
}

ACB_INLINE int
acb_contains_fmpz(const acb_t x, const fmpz_t y)
{
    return arb_contains_fmpz(acb_realref(x), y) &&
            arb_contains_zero(acb_imagref(x));
}

ACB_INLINE int
acb_contains(const acb_t x, const acb_t y)
{
    return arb_contains(acb_realref(x), acb_realref(y)) &&
            arb_contains(acb_imagref(x), acb_imagref(y));
}

ACB_INLINE int
acb_contains_interior(const acb_t x, const acb_t y)
{
    return arb_contains_interior(acb_realref(x), acb_realref(y)) &&
            arb_contains_interior(acb_imagref(x), acb_imagref(y));
}

ACB_INLINE void
acb_set_ui(acb_t z, ulong c)
{
    arb_set_ui(acb_realref(z), c);
    arb_zero(acb_imagref(z));
}

ACB_INLINE void
acb_set_d(acb_t z, double c)
{
    arb_set_d(acb_realref(z), c);
    arb_zero(acb_imagref(z));
}

ACB_INLINE void
acb_set_si(acb_t z, slong c)
{
    arb_set_si(acb_realref(z), c);
    arb_zero(acb_imagref(z));
}

ACB_INLINE void
acb_set_si_si(acb_t z, slong x, slong y)
{
    arb_set_si(acb_realref(z), x);
    arb_set_si(acb_imagref(z), y);
}

ACB_INLINE void
acb_set_d_d(acb_t z, double x, double y)
{
    arb_set_d(acb_realref(z), x);
    arb_set_d(acb_imagref(z), y);
}

ACB_INLINE void
acb_set_fmpz(acb_t z, const fmpz_t c)
{
    arb_set_fmpz(acb_realref(z), c);
    arb_zero(acb_imagref(z));
}

ACB_INLINE void
acb_set_fmpz_fmpz(acb_t z, const fmpz_t x, const fmpz_t y)
{
    arb_set_fmpz(acb_realref(z), x);
    arb_set_fmpz(acb_imagref(z), y);
}

ACB_INLINE void
acb_set_round_fmpz(acb_t z, const fmpz_t y, slong prec)
{
    arb_set_round_fmpz(acb_realref(z), y, prec);
    arb_zero(acb_imagref(z));
}

int acb_contains_int(const acb_t x);

int acb_get_unique_fmpz(fmpz_t z, const acb_t x);

ACB_INLINE void
acb_set_fmpq(acb_t z, const fmpq_t c, slong prec)
{
    arb_set_fmpq(acb_realref(z), c, prec);
    arb_zero(acb_imagref(z));
}

ACB_INLINE void
acb_set_arb(acb_t z, const arb_t c)
{
    arb_set(acb_realref(z), c);
    arb_zero(acb_imagref(z));
}

ACB_INLINE void
acb_set_arb_arb(acb_t z, const arb_t x, const arb_t y)
{
    arb_set(acb_realref(z), x);
    arb_set(acb_imagref(z), y);
}

ACB_INLINE void
acb_set_round_arb(acb_t z, const arb_t x, slong prec)
{
    arb_set_round(acb_realref(z), x, prec);
    arb_zero(acb_imagref(z));
}

ACB_INLINE void
acb_trim(acb_t z, const acb_t x)
{
    arb_trim(acb_realref(z), acb_realref(x));
    arb_trim(acb_imagref(z), acb_imagref(x));
}

ACB_INLINE void
acb_add_error_arf(acb_t x, const arf_t err)
{
    arb_add_error_arf(acb_realref(x), err);
    arb_add_error_arf(acb_imagref(x), err);
}

ACB_INLINE void
acb_add_error_mag(acb_t x, const mag_t err)
{
    arb_add_error_mag(acb_realref(x), err);
    arb_add_error_mag(acb_imagref(x), err);
}

ACB_INLINE void
acb_add_error_arb(acb_t x, const arb_t err)
{
    arb_add_error(acb_realref(x), err);
    arb_add_error(acb_imagref(x), err);
}

void acb_get_mag(mag_t z, const acb_t x);

void acb_get_mag_lower(mag_t z, const acb_t x);

void acb_get_abs_ubound_arf(arf_t u, const acb_t z, slong prec);
void acb_get_abs_lbound_arf(arf_t u, const acb_t z, slong prec);
void acb_get_rad_ubound_arf(arf_t u, const acb_t z, slong prec);

ACB_INLINE void
acb_union(acb_t res, const acb_t x, const acb_t y, slong prec)
{
    arb_union(acb_realref(res), acb_realref(x), acb_realref(y), prec);
    arb_union(acb_imagref(res), acb_imagref(x), acb_imagref(y), prec);
}

void acb_arg(arb_t r, const acb_t z, slong prec);

void acb_sgn(acb_t res, const acb_t z, slong prec);

void acb_csgn(arb_t res, const acb_t z);

void acb_real_abs(acb_t res, const acb_t z, int analytic, slong prec);
void acb_real_sgn(acb_t res, const acb_t z, int analytic, slong prec);
void acb_real_heaviside(acb_t res, const acb_t z, int analytic, slong prec);
void acb_real_floor(acb_t res, const acb_t z, int analytic, slong prec);
void acb_real_ceil(acb_t res, const acb_t z, int analytic, slong prec);
void acb_real_max(acb_t res, const acb_t x, const acb_t y, int analytic, slong prec);
void acb_real_min(acb_t res, const acb_t x, const acb_t y, int analytic, slong prec);
void acb_real_sqrtpos(acb_t res, const acb_t z, int analytic, slong prec);

void acb_sqrt_analytic(acb_t res, const acb_t z, int analytic, slong prec);
void acb_rsqrt_analytic(acb_t res, const acb_t z, int analytic, slong prec);
void acb_log_analytic(acb_t res, const acb_t z, int analytic, slong prec);
void acb_pow_analytic(acb_t res, const acb_t z, const acb_t w, int analytic, slong prec);

ACB_INLINE void
acb_add(acb_t z, const acb_t x, const acb_t y, slong prec)
{
    arb_add(acb_realref(z), acb_realref(x), acb_realref(y), prec);
    arb_add(acb_imagref(z), acb_imagref(x), acb_imagref(y), prec);
}

ACB_INLINE void
acb_sub(acb_t z, const acb_t x, const acb_t y, slong prec)
{
    arb_sub(acb_realref(z), acb_realref(x), acb_realref(y), prec);
    arb_sub(acb_imagref(z), acb_imagref(x), acb_imagref(y), prec);
}

ACB_INLINE void
acb_add_si(acb_t z, const acb_t x, slong c, slong prec)
{
    arb_add_si(acb_realref(z), acb_realref(x), c, prec);
    arb_set_round(acb_imagref(z), acb_imagref(x), prec);
}

ACB_INLINE void
acb_add_ui(acb_t z, const acb_t x, ulong c, slong prec)
{
    arb_add_ui(acb_realref(z), acb_realref(x), c, prec);
    arb_set_round(acb_imagref(z), acb_imagref(x), prec);
}

ACB_INLINE void
acb_sub_si(acb_t z, const acb_t x, slong c, slong prec)
{
    arb_sub_si(acb_realref(z), acb_realref(x), c, prec);
    arb_set_round(acb_imagref(z), acb_imagref(x), prec);
}

ACB_INLINE void
acb_sub_ui(acb_t z, const acb_t x, ulong c, slong prec)
{
    arb_sub_ui(acb_realref(z), acb_realref(x), c, prec);
    arb_set_round(acb_imagref(z), acb_imagref(x), prec);
}

ACB_INLINE void
acb_add_fmpz(acb_t z, const acb_t x, const fmpz_t y, slong prec)
{
    arb_add_fmpz(acb_realref(z), acb_realref(x), y, prec);
    arb_set_round(acb_imagref(z), acb_imagref(x), prec);
}

ACB_INLINE void
acb_add_arb(acb_t z, const acb_t x, const arb_t y, slong prec)
{
    arb_add(acb_realref(z), acb_realref(x), y, prec);
    arb_set_round(acb_imagref(z), acb_imagref(x), prec);
}

ACB_INLINE void
acb_sub_fmpz(acb_t z, const acb_t x, const fmpz_t y, slong prec)
{
    arb_sub_fmpz(acb_realref(z), acb_realref(x), y, prec);
    arb_set_round(acb_imagref(z), acb_imagref(x), prec);
}

ACB_INLINE void
acb_sub_arb(acb_t z, const acb_t x, const arb_t y, slong prec)
{
    arb_sub(acb_realref(z), acb_realref(x), y, prec);
    arb_set_round(acb_imagref(z), acb_imagref(x), prec);
}

ACB_INLINE void
acb_neg(acb_t z, const acb_t x)
{
    arb_neg(acb_realref(z), acb_realref(x));
    arb_neg(acb_imagref(z), acb_imagref(x));
}

ACB_INLINE void
acb_conj(acb_t z, const acb_t x)
{
    arb_set(acb_realref(z), acb_realref(x));
    arb_neg(acb_imagref(z), acb_imagref(x));
}

ACB_INLINE void
acb_abs(arb_t u, const acb_t z, slong prec)
{
    arb_hypot(u, acb_realref(z), acb_imagref(z), prec);
}

ACB_INLINE void
acb_mul_ui(acb_t z, const acb_t x, ulong y, slong prec)
{
    arb_mul_ui(acb_realref(z), acb_realref(x), y, prec);
    arb_mul_ui(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_mul_si(acb_t z, const acb_t x, slong y, slong prec)
{
    arb_mul_si(acb_realref(z), acb_realref(x), y, prec);
    arb_mul_si(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_mul_fmpz(acb_t z, const acb_t x, const fmpz_t y, slong prec)
{
    arb_mul_fmpz(acb_realref(z), acb_realref(x), y, prec);
    arb_mul_fmpz(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_mul_arb(acb_t z, const acb_t x, const arb_t y, slong prec)
{
    arb_mul(acb_realref(z), acb_realref(x), y, prec);
    arb_mul(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_mul_onei(acb_t z, const acb_t x)
{
    if (z == x)
    {
        arb_swap(acb_realref(z), acb_imagref(z));
        arb_neg(acb_realref(z), acb_realref(z));
    }
    else
    {
        arb_neg(acb_realref(z), acb_imagref(x));
        arb_set(acb_imagref(z), acb_realref(x));
    }
}

ACB_INLINE void
acb_div_onei(acb_t z, const acb_t x)
{
    if (z == x)
    {
        arb_swap(acb_realref(z), acb_imagref(z));
        arb_neg(acb_imagref(z), acb_imagref(z));
    }
    else
    {
        arb_set(acb_realref(z), acb_imagref(x));
        arb_neg(acb_imagref(z), acb_realref(x));
    }
}

ACB_INLINE void
acb_mul_i_pow_si(acb_t z, const acb_t x, slong k)
{
    k &= 3;

    if (k == 0)
        acb_set(z, x);
    else if (k == 1)
        acb_mul_onei(z, x);
    else if (k == 2)
        acb_neg(z, x);
    else
        acb_div_onei(z, x);
}

void acb_mul(acb_t z, const acb_t x, const acb_t y, slong prec);

void acb_mul_naive(acb_t z, const acb_t x, const acb_t y, slong prec);

ACB_INLINE void
acb_mul_2exp_si(acb_t z, const acb_t x, slong e)
{
    arb_mul_2exp_si(acb_realref(z), acb_realref(x), e);
    arb_mul_2exp_si(acb_imagref(z), acb_imagref(x), e);
}

ACB_INLINE void
acb_mul_2exp_fmpz(acb_t z, const acb_t x, const fmpz_t c)
{
    arb_mul_2exp_fmpz(acb_realref(z), acb_realref(x), c);
    arb_mul_2exp_fmpz(acb_imagref(z), acb_imagref(x), c);
}

void acb_addmul(acb_t z, const acb_t x, const acb_t y, slong prec);

void acb_submul(acb_t z, const acb_t x, const acb_t y, slong prec);

ACB_INLINE void
acb_addmul_ui(acb_t z, const acb_t x, ulong y, slong prec)
{
    arb_addmul_ui(acb_realref(z), acb_realref(x), y, prec);
    arb_addmul_ui(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_addmul_si(acb_t z, const acb_t x, slong y, slong prec)
{
    arb_addmul_si(acb_realref(z), acb_realref(x), y, prec);
    arb_addmul_si(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_submul_ui(acb_t z, const acb_t x, ulong y, slong prec)
{
    arb_submul_ui(acb_realref(z), acb_realref(x), y, prec);
    arb_submul_ui(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_submul_si(acb_t z, const acb_t x, slong y, slong prec)
{
    arb_submul_si(acb_realref(z), acb_realref(x), y, prec);
    arb_submul_si(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_addmul_fmpz(acb_t z, const acb_t x, const fmpz_t y, slong prec)
{
    arb_addmul_fmpz(acb_realref(z), acb_realref(x), y, prec);
    arb_addmul_fmpz(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_submul_fmpz(acb_t z, const acb_t x, const fmpz_t y, slong prec)
{
    arb_submul_fmpz(acb_realref(z), acb_realref(x), y, prec);
    arb_submul_fmpz(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_addmul_arb(acb_t z, const acb_t x, const arb_t y, slong prec)
{
    arb_addmul(acb_realref(z), acb_realref(x), y, prec);
    arb_addmul(acb_imagref(z), acb_imagref(x), y, prec);
}

ACB_INLINE void
acb_submul_arb(acb_t z, const acb_t x, const arb_t y, slong prec)
{
    arb_submul(acb_realref(z), acb_realref(x), y, prec);
    arb_submul(acb_imagref(z), acb_imagref(x), y, prec);
}

void acb_dot_simple(acb_t res, const acb_t initial, int subtract,
    acb_srcptr x, slong xstep, acb_srcptr y, slong ystep, slong len, slong prec);
void acb_dot_precise(acb_t res, const acb_t initial, int subtract,
    acb_srcptr x, slong xstep, acb_srcptr y, slong ystep, slong len, slong prec);
void acb_dot(acb_t res, const acb_t initial, int subtract,
    acb_srcptr x, slong xstep, acb_srcptr y, slong ystep, slong len, slong prec);

void acb_approx_dot(acb_t res, const acb_t initial, int subtract,
    acb_srcptr x, slong xstep, acb_srcptr y, slong ystep, slong len, slong prec);

void acb_dot_ui(acb_t res, const acb_t initial, int subtract,
    acb_srcptr x, slong xstep, const ulong * y, slong ystep, slong len, slong prec);
void acb_dot_si(acb_t res, const acb_t initial, int subtract,
    acb_srcptr x, slong xstep, const slong * y, slong ystep, slong len, slong prec);
void acb_dot_uiui(acb_t res, const acb_t initial, int subtract,
    acb_srcptr x, slong xstep, const ulong * y, slong ystep, slong len, slong prec);
void acb_dot_siui(acb_t res, const acb_t initial, int subtract,
    acb_srcptr x, slong xstep, const ulong * y, slong ystep, slong len, slong prec);
void acb_dot_fmpz(acb_t res, const acb_t initial, int subtract,
    acb_srcptr x, slong xstep, const fmpz * y, slong ystep, slong len, slong prec);

void acb_inv(acb_t z, const acb_t x, slong prec);

void acb_div(acb_t z, const acb_t x, const acb_t y, slong prec);

ACB_INLINE void
acb_div_ui(acb_t z, const acb_t x, ulong c, slong prec)
{
    arb_div_ui(acb_realref(z), acb_realref(x), c, prec);
    arb_div_ui(acb_imagref(z), acb_imagref(x), c, prec);
}

ACB_INLINE void
acb_div_si(acb_t z, const acb_t x, slong c, slong prec)
{
    arb_div_si(acb_realref(z), acb_realref(x), c, prec);
    arb_div_si(acb_imagref(z), acb_imagref(x), c, prec);
}

ACB_INLINE void
acb_div_arb(acb_t z, const acb_t x, const arb_t c, slong prec)
{
    arb_div(acb_realref(z), acb_realref(x), c, prec);
    arb_div(acb_imagref(z), acb_imagref(x), c, prec);
}

ACB_INLINE void
acb_div_fmpz(acb_t z, const acb_t x, const fmpz_t c, slong prec)
{
    arb_div_fmpz(acb_realref(z), acb_realref(x), c, prec);
    arb_div_fmpz(acb_imagref(z), acb_imagref(x), c, prec);
}

void acb_cube(acb_t y, const acb_t x, slong prec);
void acb_pow_fmpz(acb_t y, const acb_t b, const fmpz_t e, slong prec);
void acb_pow_ui(acb_t y, const acb_t b, ulong e, slong prec);
void acb_pow_si(acb_t y, const acb_t b, slong e, slong prec);

ACB_INLINE void
acb_const_pi(acb_t x, slong prec)
{
    arb_const_pi(acb_realref(x), prec);
    arb_zero(acb_imagref(x));
}

void acb_log(acb_t r, const acb_t z, slong prec);
void acb_log1p(acb_t r, const acb_t z, slong prec);

void acb_exp(acb_t r, const acb_t z, slong prec);
void acb_exp_pi_i(acb_t r, const acb_t z, slong prec);
void acb_exp_invexp(acb_t r, acb_t s, const acb_t z, slong prec);
void acb_expm1(acb_t r, const acb_t z, slong prec);

void acb_sin(acb_t r, const acb_t z, slong prec);
void acb_cos(acb_t r, const acb_t z, slong prec);
void acb_sin_cos(acb_t s, acb_t c, const acb_t z, slong prec);
void acb_tan(acb_t r, const acb_t z, slong prec);
void acb_cot(acb_t r, const acb_t z, slong prec);

void acb_asin(acb_t r, const acb_t z, slong prec);
void acb_acos(acb_t r, const acb_t z, slong prec);
void acb_atan(acb_t r, const acb_t z, slong prec);
void acb_asinh(acb_t r, const acb_t z, slong prec);
void acb_acosh(acb_t r, const acb_t z, slong prec);
void acb_atanh(acb_t r, const acb_t z, slong prec);

ACB_INLINE void
acb_sinh(acb_t y, const acb_t x, slong prec)
{
    acb_mul_onei(y, x);
    acb_sin(y, y, prec);
    acb_div_onei(y, y);
}

ACB_INLINE void
acb_cosh(acb_t y, const acb_t x, slong prec)
{
    acb_mul_onei(y, x);
    acb_cos(y, y, prec);
}

ACB_INLINE void
acb_sinh_cosh(acb_t y, acb_t z, const acb_t x, slong prec)
{
    acb_mul_onei(y, x);
    acb_sin_cos(y, z, y, prec);
    acb_div_onei(y, y);
}

ACB_INLINE void
acb_tanh(acb_t y, const acb_t x, slong prec)
{
    acb_mul_onei(y, x);
    acb_tan(y, y, prec);
    acb_div_onei(y, y);
}

ACB_INLINE void
acb_coth(acb_t y, const acb_t x, slong prec)
{
    acb_mul_onei(y, x);
    acb_cot(y, y, prec);
    acb_mul_onei(y, y);
}

void acb_sech(acb_t r, const acb_t z, slong prec);
void acb_csch(acb_t r, const acb_t z, slong prec);

ACB_INLINE void
acb_sec(acb_t y, const acb_t x, slong prec)
{
    acb_mul_onei(y, x);
    acb_sech(y, y, prec);
}

ACB_INLINE void
acb_csc(acb_t y, const acb_t x, slong prec)
{
    acb_mul_onei(y, x);
    acb_csch(y, y, prec);
    acb_mul_onei(y, y);
}

void acb_sin_pi(acb_t r, const acb_t z, slong prec);
void acb_cos_pi(acb_t r, const acb_t z, slong prec);
void acb_sin_cos_pi(acb_t s, acb_t c, const acb_t z, slong prec);
void acb_tan_pi(acb_t r, const acb_t z, slong prec);
void acb_cot_pi(acb_t r, const acb_t z, slong prec);
void acb_csc_pi(acb_t y, const acb_t x, slong prec);

void acb_sinc(acb_t res, const acb_t z, slong prec);
void acb_sinc_pi(acb_t res, const acb_t z, slong prec);

void acb_pow_arb(acb_t z, const acb_t x, const arb_t y, slong prec);
void acb_pow(acb_t r, const acb_t x, const acb_t y, slong prec);

void acb_sqrt(acb_t y, const acb_t x, slong prec);
void acb_rsqrt(acb_t y, const acb_t x, slong prec);
void acb_sqrts(acb_t y1, acb_t y2, const acb_t x, slong prec);

void acb_root_ui(acb_t y, const acb_t x, ulong k, slong prec);

void acb_quadratic_roots_fmpz(acb_t r1, acb_t r2,
    const fmpz_t a, const fmpz_t b, const fmpz_t c, slong prec);

void acb_chebyshev_t_ui(acb_t a, ulong n, const acb_t x, slong prec);
void acb_chebyshev_t2_ui(acb_t a, acb_t b, ulong n, const acb_t x, slong prec);
void acb_chebyshev_u_ui(acb_t a, ulong n, const acb_t x, slong prec);
void acb_chebyshev_u2_ui(acb_t a, acb_t b, ulong n, const acb_t x, slong prec);

void acb_rising_ui(acb_t z, const acb_t x, ulong n, slong prec);
void acb_rising(acb_t z, const acb_t x, const acb_t n, slong prec);
void acb_rising2_ui(acb_t u, acb_t v, const acb_t x, ulong n, slong prec);

void acb_rising_ui_get_mag(mag_t bound, const acb_t s, ulong n);

void acb_gamma(acb_t y, const acb_t x, slong prec);
void acb_rgamma(acb_t y, const acb_t x, slong prec);
void acb_lgamma(acb_t y, const acb_t x, slong prec);
void acb_log_sin_pi(acb_t res, const acb_t z, slong prec);
void acb_digamma(acb_t y, const acb_t x, slong prec);
void acb_zeta(acb_t z, const acb_t s, slong prec);
void acb_hurwitz_zeta(acb_t z, const acb_t s, const acb_t a, slong prec);
void acb_polygamma(acb_t res, const acb_t s, const acb_t z, slong prec);

void acb_bernoulli_poly_ui(acb_t res, ulong n, const acb_t x, slong prec);

void acb_log_barnes_g(acb_t res, const acb_t z, slong prec);
void acb_barnes_g(acb_t res, const acb_t z, slong prec);

void acb_polylog(acb_t w, const acb_t s, const acb_t z, slong prec);
void acb_polylog_si(acb_t w, slong s, const acb_t z, slong prec);

void acb_agm1(acb_t m, const acb_t z, slong prec);
void acb_agm1_cpx(acb_ptr m, const acb_t z, slong len, slong prec);
void acb_agm(acb_t res, const acb_t a, const acb_t b, slong prec);

#define ACB_LAMBERTW_LEFT 2
#define ACB_LAMBERTW_MIDDLE 4

void acb_lambertw_asymp(acb_t res, const acb_t z, const fmpz_t k, slong L, slong M, slong prec);
int acb_lambertw_check_branch(const acb_t w, const fmpz_t k, slong prec);
void acb_lambertw_bound_deriv(mag_t res, const acb_t z, const acb_t ez1, const fmpz_t k);
void acb_lambertw(acb_t res, const acb_t z, const fmpz_t k, int flags, slong prec);

ACB_INLINE void
acb_sqr(acb_t res, const acb_t val, slong prec)
{
    acb_mul(res, val, val, prec);
}

ACB_INLINE int
acb_is_finite(const acb_t x)
{
    return arb_is_finite(acb_realref(x)) && arb_is_finite(acb_imagref(x));
}

ACB_INLINE void
acb_indeterminate(acb_t x)
{
    arb_indeterminate(acb_realref(x));
    arb_indeterminate(acb_imagref(x));
}

ACB_INLINE acb_ptr
_acb_vec_entry_ptr(acb_ptr vec, slong i)
{
    return vec + i;
}

ACB_INLINE void
_acb_vec_zero(acb_ptr A, slong n)
{
    slong i;
    for (i = 0; i < n; i++)
        acb_zero(A + i);
}

ACB_INLINE int
_acb_vec_is_zero(acb_srcptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        if (!acb_is_zero(vec + i))
            return 0;
    return 1;
}

ACB_INLINE void
_acb_vec_set(acb_ptr res, acb_srcptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_set(res + i, vec + i);
}

ACB_INLINE void
_acb_vec_set_round(acb_ptr res, acb_srcptr vec, slong len, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_set_round(res + i, vec + i, prec);
}

ACB_INLINE void
_acb_vec_swap(acb_ptr res, acb_ptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_swap(res + i, vec + i);
}

ACB_INLINE void
_acb_vec_neg(acb_ptr res, acb_srcptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_neg(res + i, vec + i);
}

ACB_INLINE void
_acb_vec_add(acb_ptr res, acb_srcptr vec1, acb_srcptr vec2, slong len, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_add(res + i, vec1 + i, vec2 + i, prec);
}

ACB_INLINE void
_acb_vec_sub(acb_ptr res, acb_srcptr vec1, acb_srcptr vec2, slong len, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_sub(res + i, vec1 + i, vec2 + i, prec);
}

ACB_INLINE void
_acb_vec_scalar_submul(acb_ptr res, acb_srcptr vec, slong len, const acb_t c, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_submul(res + i, vec + i, c, prec);
}

ACB_INLINE void
_acb_vec_scalar_addmul(acb_ptr res, acb_srcptr vec, slong len, const acb_t c, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_addmul(res + i, vec + i, c, prec);
}

ACB_INLINE void
_acb_vec_scalar_mul(acb_ptr res, acb_srcptr vec, slong len, const acb_t c, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_mul(res + i, vec + i, c, prec);
}

ACB_INLINE void
_acb_vec_scalar_mul_ui(acb_ptr res, acb_srcptr vec, slong len, ulong c, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_mul_ui(res + i, vec + i, c, prec);
}

ACB_INLINE void
_acb_vec_scalar_mul_2exp_si(acb_ptr res, acb_srcptr vec, slong len, slong c)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_mul_2exp_si(res + i, vec + i, c);
}

ACB_INLINE void
_acb_vec_scalar_mul_onei(acb_ptr res, acb_srcptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_mul_onei(res + i, vec + i);
}

ACB_INLINE void
_acb_vec_scalar_div_ui(acb_ptr res, acb_srcptr vec, slong len, ulong c, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_div_ui(res + i, vec + i, c, prec);
}

ACB_INLINE void
_acb_vec_scalar_div(acb_ptr res, acb_srcptr vec, slong len, const acb_t c, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_div(res + i, vec + i, c, prec);
}

ACB_INLINE void
_acb_vec_scalar_mul_arb(acb_ptr res, acb_srcptr vec, slong len, const arb_t c, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_mul_arb(res + i, vec + i, c, prec);
}

ACB_INLINE void
_acb_vec_scalar_div_arb(acb_ptr res, acb_srcptr vec, slong len, const arb_t c, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
    {
        arb_div(acb_realref(res + i), acb_realref(vec + i), c, prec);
        arb_div(acb_imagref(res + i), acb_imagref(vec + i), c, prec);
    }
}

ACB_INLINE void
_acb_vec_scalar_mul_fmpz(acb_ptr res, acb_srcptr vec, slong len, const fmpz_t c, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_mul_fmpz(res + i, vec + i, c, prec);
}

ACB_INLINE void
_acb_vec_scalar_div_fmpz(acb_ptr res, acb_srcptr vec, slong len, const fmpz_t c, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_div_fmpz(res + i, vec + i, c, prec);
}

ACB_INLINE void
_acb_vec_sqr(acb_ptr res, acb_srcptr vec, slong len, slong prec)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_sqr(res + i, vec + i, prec);
}

#ifdef FLINT_HAVE_FILE
void acb_fprint(FILE * file, const acb_t x);
void acb_fprintd(FILE * file, const acb_t z, slong digits);
void acb_fprintn(FILE * fp, const acb_t z, slong digits, ulong flags);
#endif

void acb_print(const acb_t x);
void acb_printd(const acb_t z, slong digits);
void acb_printn(const acb_t x, slong digits, ulong flags);

void _acb_vec_printd(acb_srcptr vec, slong len, slong ndigits);
void _acb_vec_printn(acb_srcptr vec, slong len, slong ndigits, ulong flags);

void acb_randtest(acb_t z, flint_rand_t state, slong prec, slong mag_bits);

void acb_randtest_special(acb_t z, flint_rand_t state, slong prec, slong mag_bits);

void acb_randtest_precise(acb_t z, flint_rand_t state, slong prec, slong mag_bits);

void acb_randtest_param(acb_t z, flint_rand_t state, slong prec, slong mag_bits);

void acb_urandom(acb_t z, flint_rand_t state, slong prec);

slong acb_rel_error_bits(const acb_t x);

ACB_INLINE slong
acb_rel_accuracy_bits(const acb_t x)
{
    return -acb_rel_error_bits(x);
}

slong acb_rel_one_accuracy_bits(const acb_t x);

ACB_INLINE slong
acb_bits(const acb_t x)
{
    slong b1, b2;
    b1 = arb_bits(acb_realref(x));
    b2 = arb_bits(acb_imagref(x));
    return FLINT_MAX(b1, b2);
}

ACB_INLINE int
acb_is_real(const acb_t x)
{
    return arb_is_zero(acb_imagref(x));
}

ACB_INLINE int
_acb_vec_is_real(acb_srcptr v, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (!acb_is_real(v + i))
            return 0;
    }

    return 1;
}

ACB_INLINE int
_acb_vec_is_finite(acb_srcptr vec, slong len)
{
    return _arb_vec_is_finite((arb_srcptr) vec, 2 * len);
}

ACB_INLINE int
_acb_vec_equal(acb_srcptr vec1, acb_srcptr vec2, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (!acb_equal(vec1 + i, vec2 + i))
            return 0;
    }
    return 1;
}

ACB_INLINE int
_acb_vec_overlaps(acb_srcptr vec1, acb_srcptr vec2, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (!acb_overlaps(vec1 + i, vec2 + i))
            return 0;
    }

    return 1;
}

ACB_INLINE int
_acb_vec_contains(acb_srcptr vec1, acb_srcptr vec2, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (!acb_contains(vec1 + i, vec2 + i))
            return 0;
    }

    return 1;
}

ACB_INLINE void
_acb_vec_get_real(arb_ptr re, acb_srcptr vec, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        arb_set(re + i, acb_realref(vec + i));
    }
}

ACB_INLINE void
_acb_vec_get_imag(arb_ptr im, acb_srcptr vec, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        arb_set(im + i, acb_imagref(vec + i));
    }
}

ACB_INLINE void
_acb_vec_set_real_imag(acb_ptr vec, arb_srcptr re, arb_srcptr im, slong len)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        acb_set_arb_arb(vec + i, re + i, im + i);
    }
}

ACB_INLINE slong
_acb_vec_bits(acb_srcptr vec, slong len)
{
    return _arb_vec_bits((arb_srcptr) vec, 2 * len);
}

void _acb_vec_set_powers(acb_ptr xs, const acb_t x, slong len, slong prec);

ACB_INLINE void
_acb_vec_add_error_arf_vec(acb_ptr res, arf_srcptr err, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_add_error_arf(res + i, err + i);
}

ACB_INLINE void
_acb_vec_add_error_mag_vec(acb_ptr res, mag_srcptr err, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
    {
        mag_add(arb_radref(acb_realref(res + i)),
            arb_radref(acb_realref(res + i)), err + i);
        mag_add(arb_radref(acb_imagref(res + i)),
            arb_radref(acb_imagref(res + i)), err + i);
    }
}

ACB_INLINE void
_acb_vec_indeterminate(acb_ptr vec, slong len)
{
    _arb_vec_indeterminate((arb_ptr) vec, 2 * len);
}

ACB_INLINE void
_acb_vec_trim(acb_ptr res, acb_srcptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_trim(res + i, vec + i);
}

ACB_INLINE int
_acb_vec_get_unique_fmpz_vec(fmpz * res,  acb_srcptr vec, slong len)
{
    slong i;
    for (i = 0; i < len; i++)
        if (!acb_get_unique_fmpz(res + i, vec + i))
            return 0;
    return 1;
}

/* sort complex numbers in a nice-to-display order */
void _acb_vec_sort_pretty(acb_ptr vec, slong len);

/* roots of unity */
void acb_unit_root(acb_t res, ulong order, slong prec);
void _acb_vec_unit_roots(acb_ptr z, slong order, slong len, slong prec);

ACB_INLINE slong
acb_allocated_bytes(const acb_t x)
{
    return arb_allocated_bytes(acb_realref(x)) + arb_allocated_bytes(acb_imagref(x));
}

ACB_INLINE slong
_acb_vec_allocated_bytes(acb_srcptr vec, slong len)
{
    return _arb_vec_allocated_bytes((arb_srcptr) vec, 2 * len);
}

ACB_INLINE double
_acb_vec_estimate_allocated_bytes(slong len, slong prec)
{
    return 2 * _arb_vec_estimate_allocated_bytes(len, prec);
}

#ifdef __cplusplus
}
#endif

#endif
