/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef DFLOAT_H
#define DFLOAT_H

#include <math.h>
#include "flint.h"
#include "arb_types.h"
#include "acb_types.h"
#include "gr_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Floating-point expansions and balls built on IEEE doubles.

   A dN value (N = 1, 2, 3, 4) is an unevaluated sum
   x = d[0] + d[1] + ... + d[N-1] of doubles, roughly nonoverlapping:
   |d[i+1]| is at most about ulp(d[i]).  It carries about 53 N bits
   of precision.  The head d[0] is the double nearest to x (up to
   ties), and d[0] = 0 implies x = 0.  d1 is a plain double.

   A dNb value is a ball: a dN midpoint plus a double radius, and
   represents every real number x with |x - mid| <= rad.  The
   midpoint is always finite; rad is finite and >= 0, or +inf
   (the whole real line, which is a valid enclosure of anything).
   NaN never occurs in a valid dNb.

   The ball operations are rigorous: the output ball contains the
   image of the input balls under the exact operation, for all finite
   inputs, including where the double arithmetic underflows or
   overflows.  The error bounds are loose (constant factors are not
   optimized) but exact operations stay exact: [1 +/- 0] + [1 +/- 0]
   is [2 +/- 0], and in general the radius is zero whenever the
   inputs are exact and every floating-point operation performed on
   the midpoints was exact.

   The plain dN types are approximate arithmetic in the style of
   double: results are rounded, overflow gives inf, and inf/nan
   propagate.

   Requirements: round-to-nearest IEEE double arithmetic with a
   correctly rounded fma().  Hardware FMA (x86-64-v3 or later, or
   any 64-bit ARM) is required for performance; without it fma()
   is a slow library call but the results are still correct.  The
   error-free transformations below break if the compiler contracts
   a separate multiply and add into an fma: the module's own files
   are compiled with -ffp-contract=off, and the inline code in this
   header is written so that every product whose rounding matters is
   also consumed by an fma (which stops gcc from contracting it). */

typedef struct { double d[1]; } d1_struct;
typedef struct { double d[2]; } d2_struct;
typedef struct { double d[3]; } d3_struct;
typedef struct { double d[4]; } d4_struct;

typedef struct { double d[1]; double rad; } d1b_struct;
typedef struct { double d[2]; double rad; } d2b_struct;
typedef struct { double d[3]; double rad; } d3b_struct;
typedef struct { double d[4]; double rad; } d4b_struct;

typedef d1_struct d1_t[1];
typedef d2_struct d2_t[1];
typedef d3_struct d3_t[1];
typedef d4_struct d4_t[1];
typedef d1b_struct d1b_t[1];
typedef d2b_struct d2b_t[1];
typedef d3b_struct d3b_t[1];
typedef d4b_struct d4b_t[1];

typedef d1_struct * d1_ptr;
typedef d2_struct * d2_ptr;
typedef d3_struct * d3_ptr;
typedef d4_struct * d4_ptr;
typedef d1b_struct * d1b_ptr;
typedef d2b_struct * d2b_ptr;
typedef d3b_struct * d3b_ptr;
typedef d4b_struct * d4b_ptr;

typedef const d1_struct * d1_srcptr;
typedef const d2_struct * d2_srcptr;
typedef const d3_struct * d3_srcptr;
typedef const d4_struct * d4_srcptr;
typedef const d1b_struct * d1b_srcptr;
typedef const d2b_struct * d2b_srcptr;
typedef const d3b_struct * d3b_srcptr;
typedef const d4b_struct * d4b_srcptr;

#define DFLOAT_MAX_N 4

/* ------------------------------------------------------------------ */
/* Error-free transformations and radius helpers                       */
/* ------------------------------------------------------------------ */

/* s + e = a + b exactly (no precondition) */
#define DFLOAT_TWO_SUM(s, e, a, b) \
    do { \
        double _a = (a), _b = (b), _bb; \
        (s) = _a + _b; \
        _bb = (s) - _a; \
        (e) = (_a - ((s) - _bb)) + (_b - _bb); \
    } while (0)

/* s + e = a + b exactly, provided exponent(a) >= exponent(b) or a = 0 */
#define DFLOAT_FAST_TWO_SUM(s, e, a, b) \
    do { \
        double _a = (a), _b = (b); \
        (s) = _a + _b; \
        (e) = _b - ((s) - _a); \
    } while (0)

/* p + e = a b exactly, provided a b does not underflow
   (|a b| >= 2^-969 or a b = 0) and does not overflow */
#define DFLOAT_TWO_PROD(p, e, a, b) \
    do { \
        double _a = (a), _b = (b); \
        (p) = _a * _b; \
        (e) = fma(_a, _b, -(p)); \
    } while (0)

/* Radius arithmetic.  All radius quantities are nonnegative doubles
   (or +inf), computed with round-to-nearest, and an upper bound on
   the true quantity is restored at the end by a single multiplication
   with DFLOAT_RAD_SLACK = 1 + 2^-45, which covers up to 2^7
   nearest-rounded operations on normal numbers (each has relative
   error at most 2^-53).  The products and sums below never produce
   an inexact subnormal result: sums of doubles that land in the
   subnormal range are exact, and _dfloat_rad_mul handles products
   that would underflow explicitly.  Zero times infinity is zero
   (an exact zero times any real number). */

#define DFLOAT_RAD_SLACK 0x1.000000000008p+0            /* 1 + 2^-45 */
#define DFLOAT_RAD_UNDERFLOW_THRESHOLD 0x1p-1000
#define DFLOAT_RAD_UNDERFLOW_BUMP 0x1p-1070            /* a subnormal */

FLINT_FORCE_INLINE double _dfloat_rad_mul_slow(double a, double b, double p)
{
    if (a == 0.0 || b == 0.0)
        return 0.0;
    /* p = fl(a b) may have rounded down by up to 2^-1075 (subnormal
       range) or by a relative 2^-53 (normal); the sum below is exact
       in the subnormal range and normal otherwise */
    return p + DFLOAT_RAD_UNDERFLOW_BUMP;
}

/* upper bound (up to the final slack factor) on a b for a, b >= 0 */
FLINT_FORCE_INLINE double _dfloat_rad_mul(double a, double b)
{
    double p = a * b;
    if (FLINT_LIKELY(p >= DFLOAT_RAD_UNDERFLOW_THRESHOLD))
        return p;
    return _dfloat_rad_mul_slow(a, b, p);
}

/* upper bound (up to the slack factor) on a / b for a >= 0, b > 0 */
FLINT_FORCE_INLINE double _dfloat_rad_div(double a, double b)
{
    double q = a / b;
    if (FLINT_LIKELY(q >= DFLOAT_RAD_UNDERFLOW_THRESHOLD))
        return q;
    if (a == 0.0)
        return 0.0;
    return q + DFLOAT_RAD_UNDERFLOW_BUMP;
}

/* final upper bound from an accumulated nearest-rounded radius */
FLINT_FORCE_INLINE double _dfloat_rad_finish(double r)
{
    return r * DFLOAT_RAD_SLACK;
}

/* lower bound on a positive normal quantity computed with a few
   nearest roundings (returns 0 for subnormal inputs, which is a
   valid lower bound) */
FLINT_FORCE_INLINE double _dfloat_lower(double r)
{
    double t = r * 0x1.fffffffffff00p-1;             /* 1 - 2^-45 */
    return (r >= 0x1p-1000) ? t : 0.0;
}

/* whether a ball (n midpoint components and a radius) is the whole
   real line: a component or the radius is +-inf or NaN (x * 0 is 0 for
   finite x and NaN otherwise) */
FLINT_FORCE_INLINE int _dfloat_ball_is_whole(const double * d, int n, double rad)
{
    double s = rad * 0.0;
    int i;
    for (i = 0; i < n; i++)
        s += d[i] * 0.0;
    return !(s == 0.0);
}

/* the sum of the absolute values of an N-term expansion, an upper
   bound (up to slack) on |x| */
FLINT_FORCE_INLINE double _dfloat_abs_sum(const double * d, int n)
{
    double s = 0.0;
    int i;
    for (i = 0; i < n; i++)
        s += fabs(d[i]);
    return s;
}

/* ------------------------------------------------------------------ */
/* Per-format API                                                      */
/* ------------------------------------------------------------------ */

/* The functions below exist for X = d1, d2, d3, d4 (plain) and
   XB = d1b, d2b, d3b, d4b (balls); see the documentation.  The
   generated declarations are listed explicitly for the sake of
   editors and bindings. */

#define DFLOAT_DECLARE_PLAIN(X, T, TP, TS) \
void X##_set_d(T res, double x); \
void X##_set_dd(T res, double x0, double x1); \
void X##_zero(T res); \
void X##_one(T res); \
void X##_set(T res, TS x); \
void X##_neg(T res, TS x); \
void X##_abs(T res, TS x); \
void X##_canonicalise(T res, TS x); \
void X##_add(T res, TS x, TS y); \
void X##_sub(T res, TS x, TS y); \
void X##_mul(T res, TS x, TS y); \
void X##_sqr(T res, TS x); \
void X##_mul_d(T res, TS x, double y); \
void X##_div(T res, TS x, TS y); \
void X##_inv(T res, TS x); \
void X##_sqrt(T res, TS x); \
void X##_rsqrt(T res, TS x); \
void X##_mul_2exp_si(T res, TS x, slong e); \
void X##_renorm(T res, const double * t, slong len); \
int X##_is_zero(TS x); \
int X##_is_one(TS x); \
int X##_is_finite(TS x); \
int X##_is_nan(TS x); \
int X##_equal(TS x, TS y); \
int X##_cmp(TS x, TS y); \
int X##_sgn(TS x); \
double X##_get_d(TS x); \
void X##_get_arf(arf_t res, TS x); \
void X##_set_arf(T res, const arf_t x); \
void X##_get_arb(arb_t res, TS x); \
void X##_set_arb(T res, const arb_t x); \
void X##_set_fmpz(T res, const fmpz_t x); \
void X##_set_fmpq(T res, const fmpq_t x); \
int X##_set_str(T res, const char * s); \
char * X##_get_str(TS x, slong digits); \
void X##_randtest(T res, flint_rand_t state); \
void X##_randtest_special(T res, flint_rand_t state); \
void X##_exp(T res, TS x); \
void X##_sin(T res, TS x); \
void X##_cos(T res, TS x); \
void X##_sin_cos(T sn, T cs, TS x); \
void X##_log(T res, TS x); \
void X##_atan(T res, TS x); \
void X##_expm1(T res, TS x); \
void X##_log1p(T res, TS x); \
void X##_tan(T res, TS x); \
void X##_sinh(T res, TS x); \
void X##_cosh(T res, TS x); \
void X##_tanh(T res, TS x); \
void X##_asin(T res, TS x); \
void X##_acos(T res, TS x); \
void X##_atan2(T res, TS y, TS x); \
void X##_pow(T res, TS x, TS y); \
void X##_print(TS x); \
void _##X##_vec_add(TP res, TS x, TS y, slong len); \
void _##X##_vec_sub(TP res, TS x, TS y, slong len); \
void _##X##_vec_mul(TP res, TS x, TS y, slong len); \
void _##X##_vec_mul_scalar(TP res, TS x, slong len, TS c); \
void _##X##_vec_addmul_scalar(TP res, TS x, slong len, TS c); \
void _##X##_vec_submul_scalar(TP res, TS x, slong len, TS c); \
void _##X##_vec_dot(T res, TS initial, int subtract, TS x, TS y, slong len); \
void _##X##_vec_dot_rev(T res, TS initial, int subtract, TS x, TS y, slong len); \
void _##X##_vec_exp(TP res, TS x, slong len); \
void _##X##_vec_sin(TP res, TS x, slong len); \
void _##X##_vec_cos(TP res, TS x, slong len); \
void _##X##_vec_sin_cos(TP sres, TP cres, TS x, slong len); \
void _##X##_vec_log(TP res, TS x, slong len); \
void _##X##_vec_atan(TP res, TS x, slong len); \
void _##X##_vec_div(TP res, TS x, TS y, slong len); \
void _##X##_vec_sqrt(TP res, TS x, slong len); \
void _##X##_vec_rsqrt(TP res, TS x, slong len);

#define DFLOAT_DECLARE_BALL(X, T, TP, TS) \
void X##_set_d(T res, double x); \
void X##_set_dd(T res, double x0, double x1); \
void X##_set_d_rad(T res, double x, double rad); \
void X##_zero(T res); \
void X##_one(T res); \
void X##_indeterminate(T res); \
void X##_set(T res, TS x); \
void X##_neg(T res, TS x); \
void X##_abs(T res, TS x); \
void X##_canonicalise(T res, TS x); \
void X##_add(T res, TS x, TS y); \
void X##_sub(T res, TS x, TS y); \
void X##_mul(T res, TS x, TS y); \
void X##_sqr(T res, TS x); \
void X##_mul_d(T res, TS x, double y); \
int X##_div(T res, TS x, TS y); \
int X##_inv(T res, TS x); \
int X##_sqrt(T res, TS x); \
int X##_rsqrt(T res, TS x); \
void X##_mul_2exp_si(T res, TS x, slong e); \
void X##_add_error_d(T res, double err); \
void X##_get_mid_d(double * mid, TS x); \
int X##_is_exact(TS x); \
int X##_is_finite(TS x); \
int X##_is_zero(TS x); \
int X##_contains_zero(TS x); \
int X##_contains_d(TS x, double c); \
int X##_contains(TS x, TS y); \
int X##_overlaps(TS x, TS y); \
int X##_equal(TS x, TS y); \
int X##_cmp(int * res, TS x, TS y); \
double X##_get_d(TS x); \
void X##_get_arb(arb_t res, TS x); \
void X##_set_arb(T res, const arb_t x); \
void X##_set_arf(T res, const arf_t x); \
void X##_set_fmpz(T res, const fmpz_t x); \
void X##_set_fmpq(T res, const fmpq_t x); \
int X##_set_str(T res, const char * s); \
char * X##_get_str(TS x, slong digits); \
void X##_randtest(T res, flint_rand_t state); \
void X##_randtest_special(T res, flint_rand_t state); \
void X##_exp(T res, TS x); \
void X##_sin(T res, TS x); \
void X##_cos(T res, TS x); \
void X##_sin_cos(T sn, T cs, TS x); \
int X##_log(T res, TS x); \
void X##_atan(T res, TS x); \
void X##_expm1(T res, TS x); \
int X##_log1p(T res, TS x); \
int X##_tan(T res, TS x); \
void X##_sinh(T res, TS x); \
void X##_cosh(T res, TS x); \
void X##_tanh(T res, TS x); \
int X##_asin(T res, TS x); \
int X##_acos(T res, TS x); \
void X##_atan2(T res, TS y, TS x); \
int X##_pow(T res, TS x, TS y); \
void X##_print(TS x); \
void _##X##_vec_add(TP res, TS x, TS y, slong len); \
void _##X##_vec_sub(TP res, TS x, TS y, slong len); \
void _##X##_vec_mul(TP res, TS x, TS y, slong len); \
void _##X##_vec_mul_scalar(TP res, TS x, slong len, TS c); \
void _##X##_vec_addmul_scalar(TP res, TS x, slong len, TS c); \
void _##X##_vec_submul_scalar(TP res, TS x, slong len, TS c); \
void _##X##_vec_dot(T res, TS initial, int subtract, TS x, TS y, slong len); \
void _##X##_vec_dot_rev(T res, TS initial, int subtract, TS x, TS y, slong len); \
void _##X##_vec_exp(TP res, TS x, slong len); \
void _##X##_vec_sin(TP res, TS x, slong len); \
void _##X##_vec_cos(TP res, TS x, slong len); \
void _##X##_vec_sin_cos(TP sres, TP cres, TS x, slong len); \
int _##X##_vec_log(TP res, TS x, slong len); \
void _##X##_vec_atan(TP res, TS x, slong len); \
void _##X##_vec_div(TP res, TS x, TS y, slong len); \
int _##X##_vec_sqrt(TP res, TS x, slong len); \
int _##X##_vec_rsqrt(TP res, TS x, slong len);

DFLOAT_DECLARE_PLAIN(d1, d1_t, d1_ptr, d1_srcptr)
DFLOAT_DECLARE_PLAIN(d2, d2_t, d2_ptr, d2_srcptr)
DFLOAT_DECLARE_PLAIN(d3, d3_t, d3_ptr, d3_srcptr)
DFLOAT_DECLARE_PLAIN(d4, d4_t, d4_ptr, d4_srcptr)
DFLOAT_DECLARE_BALL(d1b, d1b_t, d1b_ptr, d1b_srcptr)
DFLOAT_DECLARE_BALL(d2b, d2b_t, d2b_ptr, d2b_srcptr)
DFLOAT_DECLARE_BALL(d3b, d3b_t, d3b_ptr, d3b_srcptr)
DFLOAT_DECLARE_BALL(d4b, d4b_t, d4b_ptr, d4b_srcptr)

/* ------------------------------------------------------------------ */
/* Complex numbers                                                     */
/* ------------------------------------------------------------------ */

/* A dNc value is a complex number re + im i with dN parts, and a dNcb
   value a complex ball in rectangular form: a dNb real part and a dNb
   imaginary part, representing every re + im i with re and im in the
   respective balls. Each part follows the conventions of the real
   balls (a nonfinite component or radius makes that part the whole
   real line). The parts are stored consecutively, so that an array of
   len complex numbers is also an array of 2 len real numbers. */

typedef struct { d1_struct re; d1_struct im; } d1c_struct;
typedef struct { d2_struct re; d2_struct im; } d2c_struct;
typedef struct { d3_struct re; d3_struct im; } d3c_struct;
typedef struct { d4_struct re; d4_struct im; } d4c_struct;
typedef struct { d1b_struct re; d1b_struct im; } d1cb_struct;
typedef struct { d2b_struct re; d2b_struct im; } d2cb_struct;
typedef struct { d3b_struct re; d3b_struct im; } d3cb_struct;
typedef struct { d4b_struct re; d4b_struct im; } d4cb_struct;

typedef d1c_struct d1c_t[1];
typedef d2c_struct d2c_t[1];
typedef d3c_struct d3c_t[1];
typedef d4c_struct d4c_t[1];
typedef d1cb_struct d1cb_t[1];
typedef d2cb_struct d2cb_t[1];
typedef d3cb_struct d3cb_t[1];
typedef d4cb_struct d4cb_t[1];

typedef d1c_struct * d1c_ptr;
typedef d2c_struct * d2c_ptr;
typedef d3c_struct * d3c_ptr;
typedef d4c_struct * d4c_ptr;
typedef d1cb_struct * d1cb_ptr;
typedef d2cb_struct * d2cb_ptr;
typedef d3cb_struct * d3cb_ptr;
typedef d4cb_struct * d4cb_ptr;

typedef const d1c_struct * d1c_srcptr;
typedef const d2c_struct * d2c_srcptr;
typedef const d3c_struct * d3c_srcptr;
typedef const d4c_struct * d4c_srcptr;
typedef const d1cb_struct * d1cb_srcptr;
typedef const d2cb_struct * d2cb_srcptr;
typedef const d3cb_struct * d3cb_srcptr;
typedef const d4cb_struct * d4cb_srcptr;

#define d1c_realref(x) (&(x)->re)
#define d1c_imagref(x) (&(x)->im)
#define d2c_realref(x) (&(x)->re)
#define d2c_imagref(x) (&(x)->im)
#define d3c_realref(x) (&(x)->re)
#define d3c_imagref(x) (&(x)->im)
#define d4c_realref(x) (&(x)->re)
#define d4c_imagref(x) (&(x)->im)
#define d1cb_realref(x) (&(x)->re)
#define d1cb_imagref(x) (&(x)->im)
#define d2cb_realref(x) (&(x)->re)
#define d2cb_imagref(x) (&(x)->im)
#define d3cb_realref(x) (&(x)->re)
#define d3cb_imagref(x) (&(x)->im)
#define d4cb_realref(x) (&(x)->re)
#define d4cb_imagref(x) (&(x)->im)

/* X: the complex prefix, R: the real prefix */
#define DFLOAT_DECLARE_COMPLEX(X, T, TP, TS, RT, RS) \
void X##_zero(T res); \
void X##_one(T res); \
void X##_onei(T res); \
void X##_set(T res, TS x); \
void X##_set_d(T res, double x); \
void X##_set_d_d(T res, double x, double y); \
void X##_set_re_im(T res, RS x, RS y); \
void X##_neg(T res, TS x); \
void X##_conj(T res, TS x); \
void X##_add(T res, TS x, TS y); \
void X##_sub(T res, TS x, TS y); \
void X##_mul(T res, TS x, TS y); \
void X##_sqr(T res, TS x); \
void X##_mul_re(T res, TS x, RS y); \
void X##_mul_onei(T res, TS x); \
void X##_mul_2exp_si(T res, TS x, slong e); \
void X##_div(T res, TS x, TS y); \
void X##_inv(T res, TS x); \
void X##_sqrt(T res, TS x); \
void X##_rsqrt(T res, TS x); \
void X##_abs(RT res, TS x); \
void X##_arg(RT res, TS x); \
int X##_is_zero(TS x); \
int X##_is_one(TS x); \
int X##_is_real(TS x); \
int X##_is_finite(TS x); \
int X##_equal(TS x, TS y); \
void X##_get_acb(acb_t res, TS x); \
void X##_set_acb(T res, const acb_t x); \
char * X##_get_str(TS x, slong digits); \
void X##_print(TS x); \
void X##_randtest(T res, flint_rand_t state); \
void X##_randtest_special(T res, flint_rand_t state); \
void X##_exp(T res, TS x); \
void X##_log(T res, TS x); \
void X##_sin(T res, TS x); \
void X##_cos(T res, TS x); \
void X##_sin_cos(T sn, T cs, TS x); \
void X##_tan(T res, TS x); \
void X##_sinh(T res, TS x); \
void X##_cosh(T res, TS x); \
void X##_tanh(T res, TS x); \
void X##_pow(T res, TS x, TS y); \
void _##X##_vec_add(TP res, TS x, TS y, slong len); \
void _##X##_vec_sub(TP res, TS x, TS y, slong len); \
void _##X##_vec_mul(TP res, TS x, TS y, slong len); \
void _##X##_vec_mul_scalar(TP res, TS x, slong len, TS c); \
void _##X##_vec_addmul_scalar(TP res, TS x, slong len, TS c); \
void _##X##_vec_submul_scalar(TP res, TS x, slong len, TS c); \
void _##X##_vec_dot(T res, TS initial, int subtract, TS x, TS y, slong len); \
void _##X##_vec_dot_rev(T res, TS initial, int subtract, TS x, TS y, slong len); \
void _##X##_vec_exp(TP res, TS x, slong len);

#define DFLOAT_DECLARE_COMPLEX_BALL(X, T, TP, TS, RT, RS) \
void X##_zero(T res); \
void X##_one(T res); \
void X##_onei(T res); \
void X##_indeterminate(T res); \
void X##_set(T res, TS x); \
void X##_set_d(T res, double x); \
void X##_set_d_d(T res, double x, double y); \
void X##_set_re_im(T res, RS x, RS y); \
void X##_neg(T res, TS x); \
void X##_conj(T res, TS x); \
void X##_add(T res, TS x, TS y); \
void X##_sub(T res, TS x, TS y); \
void X##_mul(T res, TS x, TS y); \
void X##_sqr(T res, TS x); \
void X##_mul_re(T res, TS x, RS y); \
void X##_mul_onei(T res, TS x); \
void X##_mul_2exp_si(T res, TS x, slong e); \
int X##_div(T res, TS x, TS y); \
int X##_inv(T res, TS x); \
void X##_sqrt(T res, TS x); \
int X##_rsqrt(T res, TS x); \
void X##_abs(RT res, TS x); \
void X##_arg(RT res, TS x); \
int X##_is_zero(TS x); \
int X##_is_one(TS x); \
int X##_is_real(TS x); \
int X##_is_exact(TS x); \
int X##_is_finite(TS x); \
int X##_contains_zero(TS x); \
int X##_contains(TS x, TS y); \
int X##_overlaps(TS x, TS y); \
int X##_equal(TS x, TS y); \
void X##_get_acb(acb_t res, TS x); \
void X##_set_acb(T res, const acb_t x); \
char * X##_get_str(TS x, slong digits); \
void X##_print(TS x); \
void X##_randtest(T res, flint_rand_t state); \
void X##_randtest_special(T res, flint_rand_t state); \
void X##_exp(T res, TS x); \
int X##_log(T res, TS x); \
void X##_sin(T res, TS x); \
void X##_cos(T res, TS x); \
void X##_sin_cos(T sn, T cs, TS x); \
int X##_tan(T res, TS x); \
void X##_sinh(T res, TS x); \
void X##_cosh(T res, TS x); \
int X##_tanh(T res, TS x); \
int X##_pow(T res, TS x, TS y); \
void _##X##_vec_add(TP res, TS x, TS y, slong len); \
void _##X##_vec_sub(TP res, TS x, TS y, slong len); \
void _##X##_vec_mul(TP res, TS x, TS y, slong len); \
void _##X##_vec_mul_scalar(TP res, TS x, slong len, TS c); \
void _##X##_vec_addmul_scalar(TP res, TS x, slong len, TS c); \
void _##X##_vec_submul_scalar(TP res, TS x, slong len, TS c); \
void _##X##_vec_dot(T res, TS initial, int subtract, TS x, TS y, slong len); \
void _##X##_vec_dot_rev(T res, TS initial, int subtract, TS x, TS y, slong len); \
void _##X##_vec_exp(TP res, TS x, slong len);

DFLOAT_DECLARE_COMPLEX(d1c, d1c_t, d1c_ptr, d1c_srcptr, d1_t, d1_srcptr)
DFLOAT_DECLARE_COMPLEX(d2c, d2c_t, d2c_ptr, d2c_srcptr, d2_t, d2_srcptr)
DFLOAT_DECLARE_COMPLEX(d3c, d3c_t, d3c_ptr, d3c_srcptr, d3_t, d3_srcptr)
DFLOAT_DECLARE_COMPLEX(d4c, d4c_t, d4c_ptr, d4c_srcptr, d4_t, d4_srcptr)
DFLOAT_DECLARE_COMPLEX_BALL(d1cb, d1cb_t, d1cb_ptr, d1cb_srcptr, d1b_t, d1b_srcptr)
DFLOAT_DECLARE_COMPLEX_BALL(d2cb, d2cb_t, d2cb_ptr, d2cb_srcptr, d2b_t, d2b_srcptr)
DFLOAT_DECLARE_COMPLEX_BALL(d3cb, d3cb_t, d3cb_ptr, d3cb_srcptr, d3b_t, d3b_srcptr)
DFLOAT_DECLARE_COMPLEX_BALL(d4cb, d4cb_t, d4cb_ptr, d4cb_srcptr, d4b_t, d4b_srcptr)

/* Conversions between formats (rounding to a smaller N; exact
   to a larger N). */
void dfloat_set_dfloat(double * res, int nres, const double * x, int nx);
void dfloat_ball_set_ball(double * res, int nres, const double * x, int nx);

/* ------------------------------------------------------------------ */
/* Generic (runtime N) internals, shared by the slow paths             */
/* ------------------------------------------------------------------ */

/* x = t[0] + ... + t[len-1] (any doubles, roughly sorted by decreasing
   magnitude for good results) -> nonoverlapping res[0..n-1] plus
   |dropped part| added to *err (err may be NULL) */
void _dfloat_renorm(double * res, int n, double * err, double * t, int len);

/* rigorous product, quotient and square root of balls (mid x[0..n-1],
   rad rx) for inputs of any magnitude, by scaling to unit exponent
   (scaled.c); the product is the slow path of dNb_mul, the other two
   are the only paths */
void _dfloat_ball_mul_scaled(double * res, double * rrad, int n,
    const double * x, double rx, const double * y, double ry);
void _dfloat_ball_div(double * res, double * rrad, int n,
    const double * x, double rx, const double * y, double ry);
int _dfloat_ball_rsqrt(double * res, double * rrad, int n,
    const double * x, double rx);
int _dfloat_ball_sqrt(double * res, double * rrad, int n,
    const double * x, double rx);

/* the compiled kernels, dispatched on n by the generic code */
void _dfloat_canonicalise(double * res, int n, const double * x);
void _dfloat_add_n(double * res, double * err, int n, const double * x, const double * y);
void _dfloat_sub_n(double * res, double * err, int n, const double * x, const double * y);
void _dfloat_mul_n(double * res, double * err, int n, const double * x, const double * y);
void _dfloat_sqr_n(double * res, double * err, int n, const double * x);
void _dfloat_div_approx_n(double * q, int n, const double * x, const double * y);
void _dfloat_sqrt_approx_n(double * s, int n, const double * x);
void _dfloat_rsqrt_approx_n(double * r, int n, const double * x);

#define DFLOAT_DECLARE_KERNELS(X) \
void _##X##_add_tracked(double * res, double * err, const double * x, const double * y); \
void _##X##_sub_tracked(double * res, double * err, const double * x, const double * y); \
void _##X##_mul_tracked(double * res, double * err, const double * x, const double * y); \
void _##X##_sqr_tracked(double * res, double * err, const double * x); \
void _##X##_div_approx(double * q, const double * x, const double * y); \
void _##X##_sqrt_approx(double * s, const double * x); \
void _##X##_rsqrt_approx(double * r, const double * x);
DFLOAT_DECLARE_KERNELS(d1)
DFLOAT_DECLARE_KERNELS(d2)
DFLOAT_DECLARE_KERNELS(d3)
DFLOAT_DECLARE_KERNELS(d4)

/* midpoint (n doubles) to and from arf, exactly */
void _dfloat_get_arf(arf_t res, const double * x, int n);
void _dfloat_set_arf(double * res, int n, double * err, const arf_t x);
void _dfloat_get_arb(arb_t res, const double * x, int n, double rad);
void _dfloat_set_arb(double * res, double * rad, int n, const arb_t x);
int _dfloat_part_get_arb(arb_t res, const double * x, int n, int ball);
int _dfloat_set_str(double * res, double * rad, int n, const char * s);
char * _dfloat_get_str(const double * x, int n, double rad, slong digits);
char * _dfloat_complex_get_str(const double * re, const double * im, int n, double rrad, double irad, slong digits);
void _dfloat_randtest(double * res, int n, flint_rand_t state);
void _dfloat_randtest_special(double * res, int n, flint_rand_t state);
double _dfloat_randtest_rad(flint_rand_t state);

/* tables for the exponential function (exp_tables.c): 5-term
   expansions of exp(i/256) for i = -89..89, exp(j/65536) for
   j = -128..128, log 2, and 1/k! for k = 0..12, with the tail beyond
   the 5 terms at most _dfloat_exp_tab_err */
FLINT_DLL extern const double _dfloat_exp_T1[179][5];
FLINT_DLL extern const double _dfloat_exp_T2[257][5];
FLINT_DLL extern const double _dfloat_exp_E1[179][5];
FLINT_DLL extern const double _dfloat_exp_E2[257][5];
FLINT_DLL extern const double _dfloat_exp_ln2[1][5];
FLINT_DLL extern const double _dfloat_exp_ln2_42[6];
FLINT_DLL extern const double _dfloat_exp_invfac[13][5];
/* trigonometric tables (exp_tables.c): sin and cos of i/64 (i = 0..50)
   and of j/2^13 (j = 0..64), the coefficients (-1)^(k/2)/k! of the
   sine and cosine series, and pi/2 as a 33-bit head plus 5 terms */
FLINT_DLL extern const double _dfloat_trig_A[51][2][5];
FLINT_DLL extern const double _dfloat_trig_B[65][2][5];
FLINT_DLL extern const double _dfloat_trig_A1[101][2];
FLINT_DLL extern const double _dfloat_trig_coef[15][5];
FLINT_DLL extern const double _dfloat_trig_pi2_33[6];
/* logarithm and arctangent tables: the reciprocals r (level 1: the
   nearest doubles of 1/m for the 128 bit-pattern intervals of m in
   [0.70703125, 1.4140625), see dev/gen_dfloat_exp_tables.py, with r = 1 for the
   two intervals adjacent to 1; level 2: of 1/(1 + j/2^13), j =
   -64..64, index j + 64) with L = -log(r) of the rounded r; the
   coefficients (-1)^(k+1)/k of log1p (index k); atan(i/64) for i =
   0..64, pi/2, and the coefficients (-1)^k/(2k+1) of the arctangent
   series (index k) */
FLINT_DLL extern const double _dfloat_log_R1[128];
FLINT_DLL extern const double _dfloat_log_L1[128][5];
FLINT_DLL extern const double _dfloat_log_R2[129];
FLINT_DLL extern const double _dfloat_log_L2[129][5];
FLINT_DLL extern const double _dfloat_log_coef[18][5];
FLINT_DLL extern const double _dfloat_atan_C[65][5];
FLINT_DLL extern const double _dfloat_pi_2[1][5];
FLINT_DLL extern const double _dfloat_atan_coef[16][5];
/* the constants pi, Euler's constant, Catalan's constant and log 2 */
FLINT_DLL extern const double _dfloat_const[4][5];
#define DFLOAT_CONST_PI 0
#define DFLOAT_CONST_EULER 1
#define DFLOAT_CONST_CATALAN 2
#define DFLOAT_CONST_LOG2 3
FLINT_DLL extern const double _dfloat_exp_tab_err;

/* Static relative error bounds of the exponential kernels dN_exp
   (dev/dfloat_exp_bound.py): for -746 <= x[0] <= 709.9, dN_exp(x) = exp(x) (1 + d)
   with |d| <= DFLOAT_EXP_EPS_N, up to an absolute error of 2^-1070
   for subnormal results. */
#define DFLOAT_EXP_EPS_1 0x1.04p-51
#define DFLOAT_EXP_EPS_2 0x1.c6p-101
#define DFLOAT_EXP_EPS_3 0x1.f6p-152
#define DFLOAT_EXP_EPS_4 0x1.48p-202
#define DFLOAT_EXP_EPS(n) ((n) == 1 ? DFLOAT_EXP_EPS_1 : (n) == 2 ? DFLOAT_EXP_EPS_2 : (n) == 3 ? DFLOAT_EXP_EPS_3 : DFLOAT_EXP_EPS_4)

/* the exponential function: res = exp(x) as a ball, for any finite
   midpoint and radius */
void _dfloat_ball_exp(double * res, double * rrad, int n,
    const double * x, double rx);
void _dfloat_exp(double * res, int n, const double * x);
/* the fallback of the kernels (via arb), with the same guarantee */
void _dfloat_exp_slow(double * res, int n, const double * x);

/* Static error bounds of the sine and cosine kernels dN_sin, dN_cos,
   dN_sin_cos (dev/dfloat_exp_bound.py): for |x| <= 2^20 pi/2 - 1,
   |dN_sin(x) - sin(x)| <= EPS (|x| + min(|x|, 1)) and
   |dN_cos(x) - cos(x)| <= EPS (|x| + 1) with EPS = DFLOAT_TRIG_EPS_N;
   the accuracy is absolute (relative to |x|) near the zeros. */
#define DFLOAT_TRIG_EPS_1 0x1.86p-51
#define DFLOAT_TRIG_EPS_2 0x1.66p-99
#define DFLOAT_TRIG_EPS_3 0x1.62p-149
#define DFLOAT_TRIG_EPS_4 0x1.eep-201
#define DFLOAT_TRIG_EPS(n) ((n) == 1 ? DFLOAT_TRIG_EPS_1 : (n) == 2 ? DFLOAT_TRIG_EPS_2 : (n) == 3 ? DFLOAT_TRIG_EPS_3 : DFLOAT_TRIG_EPS_4)

/* dispatchers to dN_sin, dN_cos, dN_sin_cos and the ball versions
   (want: 1 sin, 2 cos, 3 both) */
void _dfloat_sin_cos(double * sn, double * cs, int n, const double * x, int want);
void _dfloat_ball_sin_cos(double * sn, double * srad, double * cs, double * crad,
    int n, const double * x, double rx, int want);
/* the fallbacks via arb (any argument) */
void _dfloat_trig_slow(double * sn, double * cs, int n, const double * x, int want);

/* Static relative error bounds of the logarithm and arctangent kernels
   dN_log, dN_atan (dev/dfloat_exp_bound.py): for 2^-1022 <= x < 2^1022 (log) and
   any finite x (atan), |dN_f(x) - f(x)| <= EPS |f(x)| with EPS =
   DFLOAT_LOG_EPS_N, DFLOAT_ATAN_EPS_N. */
#define DFLOAT_LOG_EPS_1 0x1.a4p-50
#define DFLOAT_LOG_EPS_2 0x1.9cp-100
#define DFLOAT_LOG_EPS_3 0x1.c8p-151
#define DFLOAT_LOG_EPS_4 0x1.28p-201
#define DFLOAT_LOG_EPS(n) ((n) == 1 ? DFLOAT_LOG_EPS_1 : (n) == 2 ? DFLOAT_LOG_EPS_2 : (n) == 3 ? DFLOAT_LOG_EPS_3 : DFLOAT_LOG_EPS_4)
#define DFLOAT_ATAN_EPS_1 0x1.22p-50
#define DFLOAT_ATAN_EPS_2 0x1.5p-100
#define DFLOAT_ATAN_EPS_3 0x1.68p-151
#define DFLOAT_ATAN_EPS_4 0x1.92p-202
#define DFLOAT_ATAN_EPS(n) ((n) == 1 ? DFLOAT_ATAN_EPS_1 : (n) == 2 ? DFLOAT_ATAN_EPS_2 : (n) == 3 ? DFLOAT_ATAN_EPS_3 : DFLOAT_ATAN_EPS_4)

/* Static relative error bounds of expm1 and log1p (dev/dfloat_exp_bound.py):
   |dN_expm1(x) - expm1(x)| <= EPS |expm1(x)| for -746 <= x <= 709.9 and
   |dN_log1p(x) - log1p(x)| <= EPS |log1p(x)| for 1 + x in
   [2^-1022, 2^1022). */
#define DFLOAT_EXPM1_EPS_1 0x1.1ap-49
#define DFLOAT_EXPM1_EPS_2 0x1.aap-99
#define DFLOAT_EXPM1_EPS_3 0x1.c2p-150
#define DFLOAT_EXPM1_EPS_4 0x1.2p-200
#define DFLOAT_EXPM1_EPS(n) ((n) == 1 ? DFLOAT_EXPM1_EPS_1 : (n) == 2 ? DFLOAT_EXPM1_EPS_2 : (n) == 3 ? DFLOAT_EXPM1_EPS_3 : DFLOAT_EXPM1_EPS_4)
#define DFLOAT_LOG1P_EPS_1 0x1.a4p-50
#define DFLOAT_LOG1P_EPS_2 0x1.ecp-100
#define DFLOAT_LOG1P_EPS_3 0x1.1p-150
#define DFLOAT_LOG1P_EPS_4 0x1.64p-201
#define DFLOAT_LOG1P_EPS(n) ((n) == 1 ? DFLOAT_LOG1P_EPS_1 : (n) == 2 ? DFLOAT_LOG1P_EPS_2 : (n) == 3 ? DFLOAT_LOG1P_EPS_3 : DFLOAT_LOG1P_EPS_4)

/* The sieved power sum sum_{n <= N} n^(-s) (as acb_dirichlet_powsum_sieved
   with len = 1) over generic rings: terms in the ball ring T, phases
   in the ball ring P, and optionally the composites and block sums in
   the plain ring F (NULL: all balls), whose product and sum have
   relative errors at most dm and da; M is the table bound (M >= N/3:
   the full table of f(k), k odd <= N/3; otherwise the two-phase
   version with the table bound max(M, N^(2/3)), whose second phase
   runs on flint_get_num_threads() threads; 0: the default); returns
   a gr status. */
int gr_powsum_sieved(acb_t res, const acb_t s, ulong N, gr_ctx_t T, gr_ctx_t P, gr_ctx_t F, double dm, double da, ulong M);

/* The dfloat instantiations: dfloat_powsum_sieved chooses the
   precisions for about prec bits after the binary point (returns 0 if
   out of range: |Re s| > 2^10, |Im s| >= 1e60, too high a precision,
   or tables above DFLOAT_POWSUM_MAX_BYTES); _dfloat_powsum_sieved uses terms in dTb and phases in
   dPb, (T, P) in (1, 2), (1, 3), (2, 3), (2, 4), (3, 4), with plain
   dT composites for T <= 2; the _ball version keeps all terms balls. */
int dfloat_powsum_sieved(acb_t res, const acb_t s, ulong N, slong prec);
int _dfloat_powsum_sieved(acb_t res, const acb_t s, ulong N, int tn, int pn, ulong M);
int _dfloat_powsum_sieved_ball(acb_t res, const acb_t s, ulong N, int tn, int pn, ulong M);

/* The sums over j of Platt's multi-evaluation (as _platt_smk in
   acb_dirichlet/platt_multieval.c, for all j in [1, J]): terms in d3b
   balls with phases in dPb balls (P = pn), moments in plain
   double-double / double arithmetic with a priori error bounds;
   S_table (K*A*B entries) is overwritten; M is the bound of the table
   of terms (0: the default).  Returns 1 on success. */
int _dfloat_platt_smk(acb_ptr S_table, const fmpz * smk_points, const arb_t t0, slong A, slong B, ulong J, slong K, int pn, ulong M, slong prec);
/* the same without the factor c = exp(-i t0 log sqrt(pi)), as compact
   balls: S5[5 (k N + m) + 0..4] = Re hi, Re lo, Im hi, Im lo, radius */
int _dfloat_platt_smk_dd(double * S5, const fmpz * smk_points, const arb_t t0, slong A, slong B, ulong J, slong K, int pn, ulong M);

/* the memory limit of the tables of dfloat_powsum_sieved (less in a
   32-bit address space) */
#if FLINT_BITS == 64
#define DFLOAT_POWSUM_MAX_BYTES 4e9
#else
#define DFLOAT_POWSUM_MAX_BYTES 5e8
#endif

/* dispatchers to dN_log, dN_atan and the ball versions (the ball log
   returns GR_UNABLE for a ball not contained in (0, inf)) */
void _dfloat_log(double * res, int n, const double * x);
void _dfloat_atan(double * res, int n, const double * x);
int _dfloat_ball_log(double * res, double * rrad, int n, const double * x, double rx);
void _dfloat_ball_atan(double * res, double * rrad, int n, const double * x, double rx);
/* the fallbacks via arb (any argument; the ball log returns GR_UNABLE
   when the ball is not contained in (0, inf)) */
void _dfloat_log_slow(double * res, int n, const double * x);
void _dfloat_atan_slow(double * res, int n, const double * x);
int _dfloat_ball_log_slow(double * res, double * rrad, int n, const double * x, double rx);
void _dfloat_ball_atan_slow(double * res, double * rrad, int n, const double * x, double rx);
void _dfloat_expm1_slow(double * res, int n, const double * x);
void _dfloat_log1p_slow(double * res, int n, const double * x);
int _dfloat_ball_log1p_slow(double * res, double * rrad, int n, const double * x, double rx);
void _dfloat_expm1(double * res, int n, const double * x);
void _dfloat_log1p(double * res, int n, const double * x);
void _dfloat_ball_expm1(double * res, double * rrad, int n, const double * x, double rx);
int _dfloat_ball_log1p(double * res, double * rrad, int n, const double * x, double rx);
void _dfloat_ball_trig_slow(double * sn, double * srad, double * cs, double * crad,
    int n, const double * x, double rx, int want);

/* ------------------------------------------------------------------ */
/* Generic (runtime N) complex functions                               */
/* ------------------------------------------------------------------ */

/* error modes of the ball products and quotients: the rounding error is
   tracked exactly for exact inputs (auto) or always bounded (fast) */
#define DFLOAT_MODE_AUTO 0
#define DFLOAT_MODE_FAST 2

/* The real operations of one format (dN or dNb), through which the
   complex functions that are not performance-critical are implemented
   once for all N (complex.c): a real element is n doubles, followed by
   the radius for a ball, and a complex element two real elements. The
   operations that can fail return a gr status (always GR_SUCCESS for
   the plain formats); the mode (DFLOAT_MODE_*) is ignored by the plain
   formats. The complex product and square are the format's own. */
typedef struct
{
    int n;
    int ball;
    void (*zero)(double * res);
    void (*one)(double * res);
    void (*set_d)(double * res, double x);
    void (*set_fmpz)(double * res, const fmpz_t x);
    void (*neg)(double * res, const double * x);
    void (*abs)(double * res, const double * x);
    void (*add)(double * res, const double * x, const double * y);
    void (*sub)(double * res, const double * x, const double * y);
    void (*mul)(double * res, const double * x, const double * y, int mode);
    void (*sqr)(double * res, const double * x, int mode);
    int (*div)(double * res, const double * x, const double * y, int mode);
    int (*sqrt)(double * res, const double * x, int mode);
    void (*mul_d)(double * res, const double * x, double y);
    void (*mul_2exp_si)(double * res, const double * x, slong e);
    void (*exp)(double * res, const double * x);
    void (*expm1)(double * res, const double * x);
    int (*log)(double * res, const double * x);
    int (*log1p)(double * res, const double * x);
    void (*sin_cos)(double * s, double * c, const double * x);
    int (*tan)(double * res, const double * x);
    void (*tanh)(double * res, const double * x);
    void (*atan2)(double * res, const double * y, const double * x);
    int (*contains_zero)(const double * x);
    int (*cmp_all_d)(const double * x, double c, int sign, int strict);
    int (*get_fmpz)(fmpz_t res, const double * x);
    void (*get_arb)(arb_t res, const double * x);
    void (*set_arb)(double * res, const arb_t x);
    void (*randtest)(double * res, flint_rand_t state);
    void (*randtest_special)(double * res, flint_rand_t state);
    void (*canonicalise)(double * x);
    void (*cmul)(double * res, const double * x, const double * y, int mode);
    void (*csqr)(double * res, const double * x, int mode);
}
_dfloat_ops_struct;

typedef const _dfloat_ops_struct * _dfloat_ops_t;

FLINT_DLL extern const _dfloat_ops_struct _d1_ops, _d2_ops, _d3_ops, _d4_ops;
FLINT_DLL extern const _dfloat_ops_struct _d1b_ops, _d2b_ops, _d3b_ops, _d4b_ops;

/* the operations of the format dN (ball = 0) or dNb (ball = 1) */
_dfloat_ops_t _dfloat_ops(int n, int ball);

/* The complex functions (complex.c); res, x, y are complex elements.
   The status is that of the ball function (GR_SUCCESS, GR_DOMAIN,
   GR_UNABLE); the plain versions always compute a result and return
   GR_DOMAIN only for an exact singularity (division by zero, 1/sqrt(0)). */
int _dfloat_complex_div(_dfloat_ops_t ops, double * res, const double * x, const double * y, int mode);
int _dfloat_complex_inv(_dfloat_ops_t ops, double * res, const double * x, int mode);
void _dfloat_complex_sqrt(_dfloat_ops_t ops, double * res, const double * x, int mode);
int _dfloat_complex_rsqrt(_dfloat_ops_t ops, double * res, const double * x, int mode);
void _dfloat_complex_abs(_dfloat_ops_t ops, double * res, const double * x);
void _dfloat_complex_arg(_dfloat_ops_t ops, double * res, const double * x);
void _dfloat_complex_exp(_dfloat_ops_t ops, double * res, const double * x);
int _dfloat_complex_log(_dfloat_ops_t ops, double * res, const double * x);
void _dfloat_complex_sin_cos(_dfloat_ops_t ops, double * sn, double * cs, const double * x);
void _dfloat_complex_sinh_cosh(_dfloat_ops_t ops, double * sh, double * ch, const double * x);
int _dfloat_complex_tan(_dfloat_ops_t ops, double * res, const double * x);
int _dfloat_complex_tanh(_dfloat_ops_t ops, double * res, const double * x);
int _dfloat_complex_pow(_dfloat_ops_t ops, double * res, const double * x, const double * y);
void _dfloat_complex_get_acb(_dfloat_ops_t ops, acb_t res, const double * x);
void _dfloat_complex_set_acb(_dfloat_ops_t ops, double * res, const acb_t x);
char * _dfloat_complex_get_str_ops(_dfloat_ops_t ops, const double * x, slong digits);
void _dfloat_complex_randtest(_dfloat_ops_t ops, double * res, flint_rand_t state, int special);


/* ------------------------------------------------------------------ */
/* Generic rings                                                       */
/* ------------------------------------------------------------------ */

typedef struct
{
    int n;
    int ball;
    int strong;
    int fast;
    int cplx;
}
_dfloat_ctx_struct;

#define DFLOAT_CTX(ctx) ((_dfloat_ctx_struct *)(ctx))
#define DFLOAT_CTX_N(ctx) (DFLOAT_CTX(ctx)->n)
#define DFLOAT_CTX_BALL(ctx) (DFLOAT_CTX(ctx)->ball)
#define DFLOAT_CTX_STRONG(ctx) (DFLOAT_CTX(ctx)->strong)
#define DFLOAT_CTX_FAST(ctx) (DFLOAT_CTX(ctx)->fast)
#define DFLOAT_CTX_COMPLEX(ctx) (DFLOAT_CTX(ctx)->cplx)
#define DFLOAT_CTX_PREC(ctx) (53 * DFLOAT_CTX_N(ctx))

/* flags for gr_ctx_init_dfloat */
#define DFLOAT_BALL 1      /* rigorous balls instead of approximate numbers */
#define DFLOAT_STRONG 2    /* canonical (strongly nonoverlapping) results */
#define DFLOAT_FAST 4      /* balls: always bound the rounding errors of
                              products instead of tracking them exactly
                              for exact inputs (cheaper, loses exactness) */
#define DFLOAT_COMPLEX 8   /* complex numbers (dNc, dNcb) instead of real */

/* n = 1, 2, 3, 4; flags: a bitwise or of DFLOAT_BALL, DFLOAT_STRONG,
   DFLOAT_FAST and DFLOAT_COMPLEX (so that 0 and 1 select plain real
   numbers and real balls, respectively, with the defaults) */
int gr_ctx_init_dfloat(gr_ctx_t ctx, int n, int flags);

/* whether the module works on this target (see the documentation) */
int dfloat_is_supported(void);

#ifdef __cplusplus
}
#endif

#endif
