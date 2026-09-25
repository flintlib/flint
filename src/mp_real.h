/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef MP_REAL_H
#define MP_REAL_H

#include <stdint.h>
#include "flint.h"
#include "arb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Low-level limb-radix multiple-precision real arithmetic.

   Two levels, B = 2^FLINT_BITS throughout:

   - mp_real_t, a ball (midpoint and radius) with a limb mantissa and a
     radix-B exponent; functions mp_real_*.

   - fixed-point numbers on limb arrays: (x, n) is an unsigned n-limb
     fraction x[0], ..., x[n-1] representing sum x[i] B^(i - n), so
     0 <= x < 1 with ulp B^-n.  Outputs of size n + 1 carry a units limb
     at index n, those of size n + 2 two integral limbs.  Functions
     _mp_real_*.

   Functions on limb arrays, or with unchecked preconditions, start
   with an underscore.  The limb-level functions with an ulp-accurate
   contract take ulong * err right after their outputs and write a
   bound on |res - f(x)| in output ulps (err may be NULL); the
   documentation gives the largest value each can write. */

/* mp_real_t *****************************************************************/

/* value = (-1)^negative (d, size) B^(exp - size), the mantissa (d, size)
   with a nonzero top limb, so B^(exp - 1) <= |value| < B^exp; size = 0
   encodes 0.  The ball has radius err ulps of the bottom limb,
   B^(exp - size) (B^exp for a zero mantissa).  err < B always: a
   larger radius truncates the mantissa, a radius below one ulp pads it
   with zero limbs down to its own scale.  Exact values carry err = 0
   and no low zero limbs. */
typedef struct
{
    nn_ptr d;       /* mantissa limbs */
    slong alloc;    /* allocated limbs */
    slong size;     /* limbs in use */
    int negative;   /* sign */
    slong exp;      /* radix-B exponent */
    ulong err;      /* radius in ulps of the bottom limb */
}
mp_real_struct;

typedef mp_real_struct mp_real_t[1];
typedef mp_real_struct * mp_real_ptr;
typedef const mp_real_struct * mp_real_srcptr;

/* Exponent safety.  Exponents are not checked for overflow.  An
   operand is safe if its limb exponent lies within +-MP_REAL_EXP_MAX,
   and a bit-exponent argument (mul_2exp_si, add_error_2exp_si, ...) if
   it lies within +-FLINT_BITS MP_REAL_EXP_MAX.  One operation on safe
   operands cannot overflow an exponent computation (the bound is
   WORD_MAX / 4 in bits, since several functions form FLINT_BITS * exp);
   its result may land outside the safe range. */
#define MP_REAL_EXP_MAX (WORD_MAX / (4 * FLINT_BITS))

FLINT_FORCE_INLINE int
mp_real_exp_is_safe(slong e)
{
    return e >= -MP_REAL_EXP_MAX && e <= MP_REAL_EXP_MAX;
}

FLINT_FORCE_INLINE int
mp_real_is_safe(const mp_real_t x)
{
    return mp_real_exp_is_safe(x->exp);
}

/* Working precision.  Functions take a precision n in limbs and
   truncate mantissas to about n limbs (less when an operand is accurate
   to fewer).  At n = mp_real_prec_limbs(p) resp. mp_real_prec_bits(b),
   the basic operations (add, sub, mul, div, sqrt, rsqrt, mul_ui,
   div_ui, addmul_ui) return a relative radius at most 4 eps on exact
   operands and at most 16 eps on operands of relative radius at most
   eps, where eps = B^-p resp. 2^-b (sums of like signs and differences
   without cancellation). */
#define MP_REAL_PREC_GUARD 1

FLINT_FORCE_INLINE slong
mp_real_prec_limbs(slong prec_limbs)
{
    return prec_limbs + MP_REAL_PREC_GUARD;
}

FLINT_FORCE_INLINE slong
mp_real_prec_bits(slong prec_bits)
{
    return (prec_bits + FLINT_BITS - 1) / FLINT_BITS + MP_REAL_PREC_GUARD;
}

/* memory management and assignment */
void mp_real_init(mp_real_t x);
void mp_real_clear(mp_real_t x);
void _mp_real_grow(mp_real_t x, slong len);

FLINT_FORCE_INLINE void
mp_real_fit_length(mp_real_t x, slong len)
{
    if (x->alloc < len)
        _mp_real_grow(x, len);
}

void mp_real_swap(mp_real_t x, mp_real_t y);
void mp_real_zero(mp_real_t x);
void mp_real_set(mp_real_t res, const mp_real_t x);
void mp_real_set_ui(mp_real_t x, ulong c);
void mp_real_set_si(mp_real_t x, slong c);
/* x = (p, len) 2^e exactly (p unsigned; len may include zero top
   limbs) */
void _mp_real_set_mpn_2exp(mp_real_t x, nn_srcptr p, slong len, slong e);

FLINT_FORCE_INLINE void
mp_real_neg(mp_real_t res, const mp_real_t x)
{
    if (res != x)
        mp_real_set(res, x);
    if (res->size)
        res->negative ^= 1;
}

void mp_real_mul_2exp_si(mp_real_t res, const mp_real_t x, slong e);

/* predicates and magnitudes */
int mp_real_is_zero(const mp_real_t x);
/* e with |x| < 2^e over the whole ball */
slong mp_real_abs_bound_lt_2exp_si(const mp_real_t x);
/* e with the relative radius of x (nonzero midpoint) below 2^e */
slong mp_real_rel_radius_lt_2exp_si(const mp_real_t x);

/* radius */
void mp_real_add_error_2exp_si(mp_real_t x, slong e);
void mp_real_add_rel_error_2exp_si(mp_real_t x, slong e);
/* add v ulps at the mantissa's current anchor, resp. at B^anchor (the
   former depends on the current normalization of the mantissa) */
void _mp_real_add_error_ulps(mp_real_t x, double v);
void _mp_real_add_error_ulps_at(mp_real_t x, double v, slong anchor);

/* arithmetic (outputs may alias inputs) */
void mp_real_add(mp_real_t res, const mp_real_t a, const mp_real_t b, slong n);
void mp_real_sub(mp_real_t res, const mp_real_t a, const mp_real_t b, slong n);
void mp_real_addmul_ui(mp_real_t res, const mp_real_t x, const mp_real_t y, ulong c, slong n);
void mp_real_submul_ui(mp_real_t res, const mp_real_t x, const mp_real_t y, ulong c, slong n);
void mp_real_mul(mp_real_t res, const mp_real_t a, const mp_real_t b, slong n);
void mp_real_mul_ui(mp_real_t res, const mp_real_t a, ulong c, slong n);
void mp_real_mul_complex(mp_real_t rr, mp_real_t ri, const mp_real_t ar, const mp_real_t ai, const mp_real_t br, const mp_real_t bi, slong n);
/* res = a - b c GIVEN |a - b c| < B^E (wrong if the bound fails) */
void _mp_real_submul_bounded(mp_real_t res, const mp_real_t a, const mp_real_t b, const mp_real_t c, slong E, slong n);
void mp_real_div(mp_real_t res, const mp_real_t a, const mp_real_t b, slong n);
void mp_real_div_ui(mp_real_t res, const mp_real_t a, ulong c, slong n);

/* roots and the AGM */
void mp_real_sqrt(mp_real_t res, const mp_real_t x, slong n);
void mp_real_rsqrt(mp_real_t res, const mp_real_t x, slong n);
void mp_real_rsqrt_ui(mp_real_t res, ulong c, slong n);
void mp_real_root_ui(mp_real_t res, const mp_real_t x, ulong k, slong n);
void mp_real_rroot_ui(mp_real_t res, const mp_real_t x, ulong k, slong n);
void _mp_real_root_ui_order(mp_real_t res, const mp_real_t x, ulong k, slong n, int r, int recip);
void mp_real_agm(mp_real_t res, const mp_real_t x, const mp_real_t y, slong n);
void _mp_real_agm_order(mp_real_t res, const mp_real_t x, const mp_real_t y, slong n, int m);

/* conversions */
void mp_real_get_arb(arb_t res, const mp_real_t x);
void mp_real_print(const mp_real_t x);
/* a ball in [0, 1) as an n-limb fraction (truncating); *err is a bound
   in output ulps, UWORD_MAX meaning none below 2^FLINT_BITS - 1 */
void _mp_real_get_fixed(nn_ptr y, ulong * err, const mp_real_t x, slong n);
/* y = floor(x B^n) for a ball in [0, 1) if the radius determines it:
   returns 1, else 0 */
int _mp_real_get_fixed_floor(nn_ptr y, slong n, const mp_real_t x);

/* constants ****************************************************************/

/* c = pi/4, log 2, Euler's constant, e, log 10, Catalan's constant,
   zeta(3), zeta(5), Gamma(1/3), Gamma(1/4), 2/pi.  The ball versions
   evaluate c to about n limbs; the limb versions write floor(c B^n) -- n
   limbs for c < 1 (pi/4, log 2, Euler, Catalan, 2/pi), n + 1 with a units limb for
   the others -- and set *err = 1.  cache = 0 computes from scratch;
   cache = 1 reads a per-thread cache of floors, extending it to at
   least max(n + 5, 1.5 times the cached precision) when it is short
   (the ball version then encloses the floor with a radius of one
   ulp). */
void mp_real_const_pi4(mp_real_t res, slong n, int cache);
void mp_real_const_log2(mp_real_t res, slong n, int cache);
void mp_real_const_euler(mp_real_t res, slong n, int cache);
void mp_real_const_e(mp_real_t res, slong n, int cache);
void mp_real_const_log10(mp_real_t res, slong n, int cache);
void mp_real_const_catalan(mp_real_t res, slong n, int cache);
void mp_real_const_zeta3(mp_real_t res, slong n, int cache);
void mp_real_const_zeta5(mp_real_t res, slong n, int cache);
void mp_real_const_gamma_1_3(mp_real_t res, slong n, int cache);
void mp_real_const_gamma_1_4(mp_real_t res, slong n, int cache);
void mp_real_const_2_div_pi(mp_real_t res, slong n, int cache);

void _mp_real_const_pi4(nn_ptr res, ulong * err, slong n, int cache);
void _mp_real_const_log2(nn_ptr res, ulong * err, slong n, int cache);
void _mp_real_const_euler(nn_ptr res, ulong * err, slong n, int cache);
void _mp_real_const_e(nn_ptr res, ulong * err, slong n, int cache);
void _mp_real_const_log10(nn_ptr res, ulong * err, slong n, int cache);
void _mp_real_const_catalan(nn_ptr res, ulong * err, slong n, int cache);
void _mp_real_const_zeta3(nn_ptr res, ulong * err, slong n, int cache);
void _mp_real_const_zeta5(nn_ptr res, ulong * err, slong n, int cache);
void _mp_real_const_gamma_1_3(nn_ptr res, ulong * err, slong n, int cache);
void _mp_real_const_gamma_1_4(nn_ptr res, ulong * err, slong n, int cache);
void _mp_real_const_2_div_pi(nn_ptr res, ulong * err, slong n, int cache);

/* frees this thread's cache of constants */
void _mp_real_const_clear_cache(void);

/* a constant as an arb ball of prec bits (computed with cache = 0) */
typedef void (* mp_real_const_func)(mp_real_t, slong, int);
void _mp_real_const_arb(arb_t res, mp_real_const_func f, slong prec);

/* Euler's constant with a forced choice of logarithms (tuning) */
void _mp_real_const_euler_tune(mp_real_t res, slong n, int set);

/* series ********************************************************************/

/* atan(p/q) resp. atanh(p/q), 0 <= p < q (p/q at most about 0.99) */
void _mp_real_atan_frac_bsplit(mp_real_t res, nn_srcptr p, slong pn,
    nn_srcptr q, slong qn, slong n);
void _mp_real_atanh_frac_bsplit(mp_real_t res, nn_srcptr p, slong pn,
    nn_srcptr q, slong qn, slong n);

/* S = (coefQ + coefP sum_{k>=1} P(k)/Q(k) prod_{j<k} R(j)/Q(j)) / coefD,
   raised to power = 1 or -1, in the format of y-cruncher's
   SeriesHypergeometric; integers as signed mpn (-1)^neg (d, n) */
typedef struct
{
    nn_srcptr d;
    slong n;
    int neg;
}
mp_real_hypgeom_int_struct;

typedef struct
{
    int power;
    mp_real_hypgeom_int_struct coefP, coefQ, coefD;
    const mp_real_hypgeom_int_struct * P;
    slong Plen;
    const mp_real_hypgeom_int_struct * Q;
    slong Qlen;
    const mp_real_hypgeom_int_struct * R;
    slong Rlen;
}
mp_real_hypgeom_series_struct;

void mp_real_hypgeom_series(mp_real_t res, const mp_real_hypgeom_series_struct * s, slong n);
void mp_real_hypgeom_series_int64(mp_real_t res, int power, int64_t coefP,
    int64_t coefQ, int64_t coefD, const int64_t * P, slong Plen,
    const int64_t * Q, slong Qlen, const int64_t * R, slong Rlen, slong n);

/* elementary functions *****************************************************/

/* sin x and cos x (either output may be NULL, and may alias x) to a
   relative accuracy of about 2^-prec, or less as the radius of x
   allows: the radius enters both outputs, whose working precision is
   chosen internally.  The sine keeps its relative accuracy for small
   |x|; near other zeros the loss of significance shows in the radius.
   Arguments beyond 2^max(65536, 4 prec) give [0 +- 1]. */
void mp_real_sin_cos_bits(mp_real_t res1, mp_real_t res2, const mp_real_t x, slong prec);

/* exp x, log x and atan x to a relative accuracy of about 2^-prec, or
   less as the radius of x allows (the output may alias x); the working
   precision is chosen internally, and the results keep their relative
   accuracy near exp(x) = 1, log(x) = 0 and atan(x) = 0.  exp throws
   for |x| >= 2^(FLINT_BITS - 5), where the result's exponent would
   leave the safe range; log returns 0 (with res = 0) unless the ball
   is strictly positive, else 1. */
void mp_real_exp_bits(mp_real_t res, const mp_real_t x, slong prec);
int mp_real_log_bits(mp_real_t res, const mp_real_t x, slong prec);
void mp_real_atan_bits(mp_real_t res, const mp_real_t x, slong prec);

/* elementary functions on fixed-point numbers ******************************/

/* Rectangular-splitting series for x < 2^-32: outputs (res, n + 1)
   for exp, sin, cos, sinh, cosh (the pairs allow either output to be
   NULL), (res, n) for atan, atanh.  err <= 10 (exp), 15 (the others);
   one-sided (never above the value) for exp, sinh, cosh, atanh. */
void _mp_real_exp_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n);
void _mp_real_sin_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n);
void _mp_real_cos_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n);
void _mp_real_sin_cos_rs(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n);
void _mp_real_sinh_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n);
void _mp_real_cosh_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n);
void _mp_real_sinh_cosh_rs(nn_ptr ysinh, nn_ptr ycosh, ulong * err, nn_srcptr x, slong n);
void _mp_real_atan_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n);
void _mp_real_atanh_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n);

/* Bitwise argument reduction with tabulated log(1 + 2^-i) resp.
   atan(2^-i), for any x in [0, 1); r = 0 selects a tuned reduction
   parameter (the value of the _default_r functions), otherwise r >= 32.
   Outputs (res, n + 1) for exp, sin, cos, tan, (res, n) for log1p and
   atan.  err <= 9r + 100 (exp), 3r + 64 (log1p), 6r + 128 (sin_cos),
   8r + 256 (tan), 4r + 64 (atan), r the parameter used (the tuned
   default is at most 768). */
void _mp_real_exp_bitwise_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n, int r);
void _mp_real_log1p_bitwise_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n, int r);
void _mp_real_sin_cos_bitwise_rs(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n, int r);
void _mp_real_tan_bitwise_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n, int r);
void _mp_real_atan_bitwise_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n, int r);
int _mp_real_exp_bitwise_rs_default_r(slong n);
int _mp_real_log1p_bitwise_rs_default_r(slong n);
int _mp_real_atan_bitwise_rs_default_r(slong n);
int _mp_real_trig_bitwise_rs_default_r(slong n);

/* exp(t) resp. sin(t) and g = 1 - cos(t) for a reduced t < 2^-r
   (r >= 16; r >= 32 for alg 1, 2): outputs (y, n + 1) resp. (ysin, n),
   (yg, n), both required; alg 0 = tuned, 1 = series, 2 = sinh resp.
   sine series and a square root, 3 = one bit-burst step, 4 = full
   bit-burst.  err <= 96, computed at run time. */
void _mp_real_exp_reduced(nn_ptr y, ulong * err, nn_srcptr t, slong n, flint_bitcnt_t r, int alg);
void _mp_real_sin_cos_reduced(nn_ptr ysin, nn_ptr yg, ulong * err, nn_srcptr t, slong n, flint_bitcnt_t r, int alg);

/* exp, sin and cos on [0, 1) by halving and the reduced functions, no
   tables; outputs (y, n + 1); err <= 128.  The _r versions force the
   halving depth r (tuning). */
void _mp_real_exp_notab(nn_ptr y, ulong * err, nn_srcptr x, slong n);
void _mp_real_exp_notab_r(nn_ptr y, ulong * err, nn_srcptr x, slong n, int r);
void _mp_real_sin_cos_notab(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n);
void _mp_real_sin_cos_notab_r(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n, int r);

/* exp, sin and cos, tan on [0, 1) by diophantine (multi-prime)
   argument reduction; outputs (y, n + 1); err <= 3 (exp), 4 (sin,
   cos, tan).  The _tune versions take the number of (Gaussian) primes
   and the weight budget of the reduction. */
void _mp_real_exp_diophantine(nn_ptr y, ulong * err, nn_srcptr x, slong n);
void _mp_real_exp_diophantine_tune(nn_ptr y, ulong * err, nn_srcptr x, slong n, slong num_primes, double max_weight);
void _mp_real_sin_cos_diophantine(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n);
void _mp_real_sin_cos_diophantine_tune(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n, slong num_primes, double max_weight);
void _mp_real_tan_diophantine(nn_ptr res, ulong * err, nn_srcptr x, slong n);
void _mp_real_tan_diophantine_tune(nn_ptr res, ulong * err, nn_srcptr x, slong n, slong num_primes, double max_weight);

/* -log(x) for x in [1/2, 1) and atan(x) for x in [0, 1), outputs
   (y, n): Newton-Taylor steps over the forward functions, resp. the
   Sasaki-Kanada AGM formula.  err <= 2, computed at run time.  The
   _tune versions take the forward algorithm (0 = diophantine,
   1 = bitwise, 2 = notab) and the number of series terms (0 = default)
   resp. theta terms. */
void _mp_real_neglog_newton(nn_ptr y, ulong * err, nn_srcptr x, slong n);
void _mp_real_neglog_newton_tune(nn_ptr y, ulong * err, nn_srcptr x, slong n, int forward, slong N);
void _mp_real_atan_newton(nn_ptr y, ulong * err, nn_srcptr x, slong n);
void _mp_real_atan_newton_tune(nn_ptr y, ulong * err, nn_srcptr x, slong n, int forward, slong N);
void _mp_real_neglog_agm(nn_ptr y, ulong * err, nn_srcptr x, slong n);
void _mp_real_neglog_agm_tune(nn_ptr y, ulong * err, nn_srcptr x, slong n, slong N);

/* Newton inverses and roots on fixed-point numbers *************************/

/* Not ulp-accurate; relative bounds (with (a, an) the input):
   inv: a in [1/B, 1), a_{an-1} != 0 -> (q, n + 2), |error| <= 4 B^-n / a;
   div: b in [0, 1) over a as for inv -> (q, n + 2), 4 B^-n / a;
   rsqrt_ui: 2 <= a < B -> (res, n) fraction limbs, 2 B^-n;
   rsqrt: a in [B^-2, 1), one of the top two limbs nonzero -> (q, n + 2),
   4 B^-n / sqrt(a); sqrt: as rsqrt -> (q, n + 2), 4 B^-n / sqrt(a). */
void _mp_real_inv_newton_basecase(nn_ptr q, nn_srcptr a, slong an, slong n);
void _mp_real_inv_newton(nn_ptr q, nn_srcptr a, slong an, slong n);
void _mp_real_div_newton_invmul(nn_ptr q, nn_srcptr b, slong bn, nn_srcptr a, slong an, slong n);
void _mp_real_div_newton(nn_ptr q, nn_srcptr b, slong bn, nn_srcptr a, slong an, slong n);
void _mp_real_rsqrt_ui_newton_basecase(nn_ptr res, ulong a, slong n);
void _mp_real_rsqrt_ui_newton(nn_ptr res, ulong a, slong n);
void _mp_real_rsqrt_newton_basecase(nn_ptr q, nn_srcptr a, slong an, slong n);
void _mp_real_rsqrt_newton(nn_ptr q, nn_srcptr a, slong an, slong n);
void _mp_real_sqrt_newton_rsqrtmul(nn_ptr q, nn_srcptr a, slong an, slong n);
void _mp_real_sqrt_newton(nn_ptr q, nn_srcptr a, slong an, slong n);

/* reduction tables (also used by arb) **************************************/

/* Integer relations sum_j d_ij alpha_j = epsilon_i among the logarithms
   of the first num primes (gaussian = 0) or pi/2 and the arguments
   2 arg(pi_j) of the first num nonreal Gaussian primes (gaussian = 1). */
#define MP_REAL_REL_MAX 64
#define MP_REAL_LOG_PRIMES_MAX MP_REAL_REL_MAX
#define MP_REAL_REL_TERMINATOR -32768

typedef struct
{
    slong num;
    slong rows;
    int gaussian;
    int is_static;
    const ulong * primes;        /* the rational primes; NULL if gaussian */
    const float * weights;       /* log2 p_j, or log N(pi_j); weights[0] = 0 */
    const short * d;             /* rows x num, then a row starting with
                                    MP_REAL_REL_TERMINATOR */
    const double * epsilon;
    const double * epsilon_inv;
    double epsilon_min;          /* |epsilon| of the last row */
}
mp_real_rel_struct;

const mp_real_rel_struct * _mp_real_rel_table(int gaussian, slong num);
int _mp_real_rel_table_is_cached(int gaussian, slong num);

/* Machin-type sets: log p_i resp. arg pi_i = (1/den) sum_j c_ij
   atanh(1/x_j) resp. atan(1/x_j) */
#define MP_REAL_MACHIN_X_LIMBS (128 / FLINT_BITS)
#define MP_REAL_MACHIN_MAX_NUM 64
#define MP_REAL_MACHIN_MAX_PRIMES 8

typedef struct
{
    slong num;
    const ulong * x;             /* MP_REAL_MACHIN_X_LIMBS limbs each */
    ulong den;
    int cbits;                   /* one more than the largest coefficient
                                    bit length */
    int gaussian;
}
mp_real_machin_struct;

const mp_real_machin_struct * _mp_real_machin_table(int gaussian, slong num);
slong _mp_real_machin_table_max(int gaussian);
void _mp_real_machin_get_x(fmpz_t q, const mp_real_machin_struct * tab, slong i);
void _mp_real_machin_get_c(fmpz_t c, const mp_real_machin_struct * tab, slong i, slong j);
void _mp_real_machin_get_c_row(fmpz * row, const mp_real_machin_struct * tab, slong i);
nn_srcptr _mp_real_machin_c_row_raw(const mp_real_machin_struct * tab, slong i,
    slong * climbs, const unsigned char ** csign);

/* real and imaginary parts, consecutively, of the first 64 nonreal
   Gaussian primes in order of norm: 1+i, 1+2i, 2+3i, ... (FLINT_DLL:
   the build exports functions automatically on Windows, not data) */
#define MP_REAL_ATAN_GAUSS_MAX 64
FLINT_DLL extern const signed char _mp_real_gaussian_primes[2 * MP_REAL_ATAN_GAUSS_MAX];

/* the angles pi/2, 2 arg(pi_j), j < num, as arb balls */
void _mp_real_atan_gauss_vec_arb(arb_ptr res, slong num, slong prec);

#ifdef __cplusplus
}
#endif

#endif
