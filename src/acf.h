/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACF_H
#define ACF_H

#ifdef ACF_INLINES_C
#define ACF_INLINE
#else
#define ACF_INLINE static inline
#endif

#include "arf.h"
#include "acf_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define acf_realref(x) (&(x)->real)
#define acf_imagref(x) (&(x)->imag)

ACF_INLINE void
acf_init(acf_t x)
{
    arf_init(acf_realref(x));
    arf_init(acf_imagref(x));
}

ACF_INLINE void
acf_clear(acf_t x)
{
    arf_clear(acf_realref(x));
    arf_clear(acf_imagref(x));
}

ACF_INLINE acf_ptr _acf_vec_init(slong n) { return (acf_ptr) _arf_vec_init(2 * n); }
ACF_INLINE void _acf_vec_clear(acf_ptr v, slong n) { _arf_vec_clear((arf_ptr) v, 2 * n); }

ACF_INLINE arf_ptr acf_real_ptr(acf_t z) { return acf_realref(z); }
ACF_INLINE arf_ptr acf_imag_ptr(acf_t z) { return acf_imagref(z); }

ACF_INLINE void
acf_set(acf_t z, const acf_t x)
{
    arf_set(acf_realref(z), acf_realref(x));
    arf_set(acf_imagref(z), acf_imagref(x));
}

ACF_INLINE void
acf_swap(acf_t z, acf_t x)
{
    FLINT_SWAP(acf_struct, *z, *x);
}

ACF_INLINE int
acf_equal(const acf_t x, const acf_t y)
{
    return arf_equal(acf_realref(x), acf_realref(y)) &&
           arf_equal(acf_imagref(x), acf_imagref(y));
}

/* todo: document */
ACF_INLINE void
acf_printd(const acf_t x, slong n)
{
    arf_printd(acf_realref(x), n);
    flint_printf(" + ");
    arf_printd(acf_imagref(x), n);
    flint_printf("*I");
}

/* todo: document */
ACF_INLINE slong
acf_bits(const acf_t x)
{
    slong b1, b2;
    b1 = arf_bits(acf_realref(x));
    b2 = arf_bits(acf_imagref(x));
    return FLINT_MAX(b1, b2);
}

ACF_INLINE slong
acf_allocated_bytes(const acf_t x)
{
    return arf_allocated_bytes(acf_realref(x)) + arf_allocated_bytes(acf_imagref(x));
}

/* todo: document */
ACF_INLINE void acf_randtest(acf_t x, flint_rand_t state, slong bits, slong mag_bits)
{
    arf_randtest(acf_realref(x), state, bits, mag_bits);
    arf_randtest(acf_imagref(x), state, bits, mag_bits);
}

/* todo: document */
ACF_INLINE void
acf_get_mag(mag_t res, const acf_t x)
{
    mag_t t, u;
    mag_init(t);
    mag_init(u);
    arf_get_mag(t, acf_realref(x));
    arf_get_mag(u, acf_imagref(x));
    mag_hypot(res, t, u);
    mag_clear(t);
    mag_clear(u);
}

/* todo: document */
ACF_INLINE void
acf_neg(acf_t z, const acf_t x)
{
    arf_neg(acf_realref(z), acf_realref(x));
    arf_neg(acf_imagref(z), acf_imagref(x));
}

/* todo: document */
ACF_INLINE int
acf_set_round(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
{
    return arf_set_round(acf_realref(res), acf_realref(x), prec, rnd) |
           (arf_set_round(acf_imagref(res), acf_imagref(x), prec, rnd) << 1);
}

/* todo: document */
ACF_INLINE int
acf_neg_round(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd)
{
    return arf_neg_round(acf_realref(res), acf_realref(x), prec, rnd) |
           (arf_neg_round(acf_imagref(res), acf_imagref(x), prec, rnd) << 1);
}

ACF_INLINE int
acf_add(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
{
    return arf_add(acf_realref(res), acf_realref(x), acf_realref(y), prec, rnd) |
           (arf_add(acf_imagref(res), acf_imagref(x), acf_imagref(y), prec, rnd) << 1);
}

ACF_INLINE int
acf_sub(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
{
    return arf_sub(acf_realref(res), acf_realref(x), acf_realref(y), prec, rnd) |
           (arf_sub(acf_imagref(res), acf_imagref(x), acf_imagref(y), prec, rnd) << 1);
}

ACF_INLINE int
acf_mul(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd)
{
    if (x == y)
        return arf_complex_sqr(acf_realref(res), acf_imagref(res),
            acf_realref(x), acf_imagref(x),
            prec, rnd);
    else
        return arf_complex_mul(acf_realref(res), acf_imagref(res),
            acf_realref(x), acf_imagref(x),
            acf_realref(y), acf_imagref(y),
            prec, rnd);
}

void acf_approx_inv(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
void acf_approx_div(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd);
void acf_approx_sqrt(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
void acf_approx_rsqrt(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);

/* the reciprocal, quotient and square roots with each part rounded in
   the mode rnd as for the elementary functions below (acf/inv.c,
   acf/sqrt.c): the
   reciprocal and quotient from the exact numerators and denominators,
   the roots by retries of the kernels in the strict modes; the return
   value carries the inexact flags of the parts.  A zero divisor gives
   the results of arf_div; sqrt(0) = 0 and rsqrt(0) = +inf. */
int acf_inv(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_div(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd);
int acf_sqrt(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_rsqrt(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);

/* elementary functions (acf/elem.c): in the strict modes each part
   correctly rounded in the mode rnd as the real functions of arf round
   their results; in ARF_RND_FAST resp. ARF_RND_ACCURATE each part
   within 1 resp. 0.51 ulp of the larger part of the value (the relative
   accuracy of the complex number: a much smaller part may be inaccurate
   or zero); the return value carries the inexact flags of the
   parts (bit 0 the real part, bit 1 the imaginary part; bits 2 and 3
   those of res2) in the strict modes.  A nan or infinite part gives
   nan parts; the arguments beyond the ranges of the real functions
   (exp beyond 2^max(128, 2 prec), the trigonometric functions beyond
   2^max(65536, 4 prec)) likewise.  log(0) = -inf, pow(0, y) = 0 for
   Re(y) > 0 and x^0 = 1. */
int acf_exp(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_exp_pi_i(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_log(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_sin(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_cos(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_sin_cos(acf_t res1, acf_t res2, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_sin_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_cos_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_sin_cos_pi(acf_t res1, acf_t res2, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_sinh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_cosh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_sinh_cosh(acf_t res1, acf_t res2, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_pow(acf_t res, const acf_t x, const acf_t y, slong prec, arf_rnd_t rnd);

/* tan, cot, sec, csc, their pi versions and the hyperbolic ones
   (acf/elem.c): the parts as quotients by |cos z|^2 or |sin z|^2 */
int acf_tan(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_cot(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_sec(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_csc(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_tan_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_cot_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_sec_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_csc_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_tanh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_coth(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_sech(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_csch(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);

/* the inverse functions (principal branches, the reciprocal ones as
   f(1/z) as in gr), expm1, log1p and sinc (acf/elem_acb.c): from balls of
   acb at increasing precisions, with the values on the axes whose parts
   are exactly zero or a multiple of pi/2 set directly */
int acf_expm1(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_log1p(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_sinc(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_sinc_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_asin(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_acos(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_atan(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_acot(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_asec(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_acsc(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_asin_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_acos_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_atan_pi(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_asinh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_acosh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_atanh(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_acoth(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_asech(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);
int acf_acsch(acf_t res, const acf_t x, slong prec, arf_rnd_t rnd);

void acf_approx_dot(acf_t res, const acf_t initial, int subtract, acf_srcptr x, slong xstep, acf_srcptr y, slong ystep, slong len, slong prec, arf_rnd_t rnd);

#ifdef __cplusplus
}
#endif

#endif
