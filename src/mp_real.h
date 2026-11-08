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

   - floating-point mantissas, the representation of arf and nfloat: a
     normalised limb array (d, n) (top bit set) with a radix-2 exponent
     e and a sign, standing for (-1)^sgnbit (d, n) 2^(e - n FLINT_BITS).
     Functions taking these end in _fp (_mp_real_exp_fp, ...).  They
     frame the argument themselves and return an unnormalised
     fixed-point result (y, nk + 1) at nk kernel limbs with an error
     count, plus a radix-2 shift as the return value where the result
     carries a binary exponent; the callers (nfloat, arf, arb, ...)
     round and normalise.  Temporaries are sized by nk, so any
     precision is allowed.

   Radix-B exponents are limb counts and radix-2 exponents bit counts,
   both slong.  Functions on limb arrays, or with unchecked
   preconditions, start with an underscore.  The limb-level functions
   with an ulp-accurate contract take ulong * err right after their
   outputs and write a bound on |res - f(x)| in output ulps (err may be
   NULL), two-sided unless the function documents a side; the
   documentation gives the largest value each can write.  The _fp
   functions take the absolute accuracy wanted in bits, pb = p +
   _mp_real_elem_guard_bits(p) (+ extra) for a relative accuracy of
   about 2^-p, and run their kernels at nk = _mp_real_elem_limbs(pb)
   limbs, plus the leading zero bits of a result that needs relative
   accuracy.  The guard bits cover the kernels' error bounds and the
   argument reduction: the bounds stay below 2^12 units below
   MP_REAL_ELEM_GUARD_LARGE_PREC bits and grow with the bitwise
   parameter beyond, to about 2^13 for exp and 2^14 for the hyperbolic
   functions, which take MP_REAL_ELEM_GUARD_BITS_LARGE guard bits. */

#define MP_REAL_ELEM_GUARD_BITS 12
#define MP_REAL_ELEM_GUARD_BITS_LARGE 16
#define MP_REAL_ELEM_GUARD_LARGE_PREC 4096

FLINT_FORCE_INLINE slong
_mp_real_elem_guard_bits(slong p)
{
    return (p < MP_REAL_ELEM_GUARD_LARGE_PREC) ? MP_REAL_ELEM_GUARD_BITS : MP_REAL_ELEM_GUARD_BITS_LARGE;
}

/* the kernel limbs for pb bits (at least two limbs on 32-bit machines,
   as some kernels require) */
FLINT_FORCE_INLINE slong
_mp_real_elem_limbs(slong pb)
{
    slong n = (pb + FLINT_BITS - 1) / FLINT_BITS;
#if FLINT_BITS == 32
    n = FLINT_MAX(n, 2);
#endif
    return n;
}

/* The output limbs _mp_real_log_fp, _mp_real_log_fp_huge (for a bignum
   exponent of cn limbs) and _mp_real_log1p_fp need for an argument of n
   limbs (and exponent e for log1p) at pb bits; the callers size their
   buffers with these */
FLINT_FORCE_INLINE slong
_mp_real_log_fp_ylen(slong pb, slong n)
{
    return _mp_real_elem_limbs(pb) + n + 5;
}

FLINT_FORCE_INLINE slong
_mp_real_log_fp_huge_ylen(slong pb, slong cn)
{
    return _mp_real_elem_limbs(pb) + cn + 3;
}

FLINT_FORCE_INLINE slong
_mp_real_log1p_fp_ylen(slong pb, slong n, slong e)
{
    slong nk = _mp_real_elem_limbs(pb);
    return nk + n + 6 + ((e < 0) ? -e : FLINT_MIN(e, FLINT_BITS * (nk + 3))) / FLINT_BITS;
}

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
void mp_real_set_fmpz(mp_real_t x, const fmpz_t f);
void mp_real_set_trunc(mp_real_t res, const mp_real_t x, slong n);
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
/* res = x^e by binary powering (x^0 = 1) */
void mp_real_pow_ui(mp_real_t res, const mp_real_t x, ulong e, slong n);

/* roots and the AGM */
void mp_real_sqrt(mp_real_t res, const mp_real_t x, slong n);
void mp_real_rsqrt(mp_real_t res, const mp_real_t x, slong n);
void mp_real_rsqrt_ui(mp_real_t res, ulong c, slong n);
/* k-th roots for 1 <= k < MP_REAL_ROOT_K_MAX (the exact series
   coefficients use the words bk + i k, i < 16, and k fits a slong) */
#if FLINT_BITS == 64
#define MP_REAL_ROOT_K_MAX (UWORD(1) << 40)
#else
#define MP_REAL_ROOT_K_MAX (UWORD(1) << 27)
#endif

void mp_real_root_ui(mp_real_t res, const mp_real_t x, ulong k, slong n);
void mp_real_rroot_ui(mp_real_t res, const mp_real_t x, ulong k, slong n);
void _mp_real_root_ui_order(mp_real_t res, const mp_real_t x, ulong k, slong n, int r, int recip);
void mp_real_agm(mp_real_t res, const mp_real_t x, const mp_real_t y, slong n);
void _mp_real_agm_order(mp_real_t res, const mp_real_t x, const mp_real_t y, slong n, int m);

/* conversions */
void mp_real_get_arb(arb_t res, const mp_real_t x);
/* the unique integer in the ball x (exactly, as an mp_real, without the
   size limits of fmpz) and 1, or 0 (res unchanged) if there is none or
   more than one; res may alias x (in place, no reallocation) */
int mp_real_unique_integer(mp_real_t res, const mp_real_t x);

/* p(n) exactly, n = nhi B + nlo (the partition function: lookup table,
   pentagonal recurrence, then the Hardy-Ramanujan-Rademacher formula);
   _mp_real_partitions_hrr by the formula for n >= 2 (p(n) for any fmpz n), with flags
   MP_REAL_PARTITIONS_NO_DFLOAT (no dfloat batches) and
   MP_REAL_PARTITIONS_TEST (the rarely taken paths at any n, for tests) */
#define MP_REAL_PARTITIONS_NO_DFLOAT 1
#define MP_REAL_PARTITIONS_TEST 2
void mp_real_partitions_hrr(mp_real_t res, ulong nhi, ulong nlo);
void _mp_real_partitions_hrr(mp_real_t res, const fmpz_t n, int flags);
/* dfloat expansions as double arrays: res = sum x[i] +- rad (exact up to
   the radius's rounding; throws for nonfinite input); and the greedy
   rounding of x to n <= 8 doubles with a radius (rad may be NULL), the
   whole line [0 +- inf] beyond the double range */
void mp_real_set_dfloat(mp_real_t res, const double * x, slong n, double rad);
void mp_real_get_dfloat(double * res, double * rad, slong n, const mp_real_t x);
void mp_real_print(const mp_real_t x);
/* a ball in [0, 1) as an n-limb fraction (truncating); *err is a bound
   in output ulps, UWORD_MAX meaning none below 2^FLINT_BITS - 1 */
void _mp_real_get_fixed(nn_ptr y, ulong * err, const mp_real_t x, slong n);
/* y = floor(x B^n) for a ball in [0, 1) if the radius determines it:
   returns 1, else 0 */
int _mp_real_get_fixed_floor(nn_ptr y, slong n, const mp_real_t x);

/* constants ****************************************************************/

/* c = pi/4, log 2, Euler's constant, e, log 10, Catalan's constant,
   zeta(3), zeta(5), Gamma(1/3), Gamma(1/4), 2/pi, 1/log 2 - 1,
   Khinchin's constant.  The ball versions evaluate c to about n limbs;
   the limb versions write floor(c B^n) -- n limbs for c < 1 (pi/4,
   log 2, Euler, Catalan, 2/pi, 1/log 2 - 1), n + 1 with a units limb
   for the others -- and set *err = 1.  cache = 0 computes from scratch;
   cache = 1 reads a shared cache of floors (one per constant, read by
   every thread; thread_support.h), extending it to at least
   max(n + 5, 1.5 times the cached precision) when it is short
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
void mp_real_const_inv_log2_frac(mp_real_t res, slong n, int cache);
void mp_real_const_khinchin(mp_real_t res, slong n, int cache);

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
void _mp_real_const_inv_log2_frac(nn_ptr res, ulong * err, slong n, int cache);
void _mp_real_const_khinchin(nn_ptr res, ulong * err, slong n, int cache);

/* frees the shared caches of the constants (no other thread may be
   using them; flint_cleanup_master does the same) */
void _mp_real_const_clear_cache(void);

/* verbose progress reports on stderr for huge computations: 0 silent
   (the default, unless the environment variable FLINT_MP_REAL_VERBOSE
   sets another level), 1 a line every few seconds from the long
   loops, 2 every phase with its size and time */
void mp_real_set_verbose(int level);
int mp_real_get_verbose(void);

/* a constant as an arb ball of prec bits (computed with cache = 0) */
typedef void (* mp_real_const_func)(mp_real_t, slong, int);
void _mp_real_const_arb(arb_t res, mp_real_const_func f, slong prec);

/* Euler's constant with a forced choice of logarithms (tuning) */
void _mp_real_const_euler_tune(mp_real_t res, slong n, int set);

/* Khinchin's constant with forced parameters K0 (logarithm rows below
   K0), n1 (Bernoulli columns below n1) and nb (Bernoulli chunks run at
   once); 0 = default.  They are clamped to 3 <= K0 <= 300,
   max(n + 1, FLINT_BITS n / 40 + 1) <= n1 <= M + 1 (M the last column
   for K0) and
   1 <= nb <= the number of threads */
void _mp_real_const_khinchin_tune(mp_real_t res, slong n, slong K0, slong n1, slong nb);
/* the default K0, n1 and nb at n limbs (nb for the current number of
   threads) */
void _mp_real_const_khinchin_params(slong * K0, slong * n1, slong * nb, slong n);

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

/* tan x to a relative accuracy of about 2^-prec, or less as the radius
   of x allows (the output may alias x); keeps its relative accuracy
   near the zeros and poles of the tangent.  Returns 1, or 0 with res
   set to zero when the ball may contain a pole (or its argument is
   beyond 2^max(65536, 4 prec)). */
int mp_real_tan_bits(mp_real_t res, const mp_real_t x, slong prec);

/* sin(pi x), cos(pi x) and tan(pi x) as the functions above, with an
   exact argument reduction: the outputs keep their relative accuracy
   near every zero (and tan near every pole), at any magnitude of x,
   and the values at exact multiples of 1/2 are exact (0, +-1), as is
   tan at the odd multiples of 1/4 (+-1). */
void mp_real_sin_cos_pi_bits(mp_real_t res1, mp_real_t res2, const mp_real_t x, slong prec);
int mp_real_tan_pi_bits(mp_real_t res, const mp_real_t x, slong prec);

/* sin(pi p/q) and cos(pi p/q) for words p, q >= 1 (either output may be
   NULL; they must be distinct) to a relative accuracy of about B^-n;
   the rational values 0, +-1, +-1/2 are exact.  tan(pi p/q) likewise,
   returning 0 (with res set to zero) at a pole, else 1; 0, +-1 are
   exact. */
void mp_real_sin_cos_pi_ui_div_ui(mp_real_t res1, mp_real_t res2, ulong p, ulong q, slong n);
int mp_real_tan_pi_ui_div_ui(mp_real_t res, ulong p, ulong q, slong n);

/* exp x, log x and atan x to a relative accuracy of about 2^-prec, or
   less as the radius of x allows (the output may alias x); the working
   precision is chosen internally, and the results keep their relative
   accuracy near exp(x) = 1, log(x) = 0 and atan(x) = 0.  exp throws
   for |x| >= 2^(FLINT_BITS - 5), where the result's exponent would
   leave the safe range; log returns 0 (with res = 0) unless the ball
   is strictly positive, else 1. */
void mp_real_exp_bits(mp_real_t res, const mp_real_t x, slong prec);
void mp_real_exp_notab_log2(mp_real_t res, const mp_real_t x, slong n);
void mp_real_exp_notab_squaring(mp_real_t res, const mp_real_t x, slong n);
void mp_real_exp_agm(mp_real_t res, const mp_real_t x, slong n);
int mp_real_log_bits(mp_real_t res, const mp_real_t x, slong prec);
void mp_real_atan_bits(mp_real_t res, const mp_real_t x, slong prec);

/* log1p x, expm1 x, sinh x and cosh x (either output may be NULL),
   tanh x, asin x, acos x, asinh x, acosh x and atanh x to a relative
   accuracy of about 2^-prec, or less as the radius of x allows (the
   outputs may alias x), as compositions of exp, log and atan without
   cancellation: the results keep their relative accuracy near their
   zeros, and small and large arguments are handled without squaring
   them.  log1p returns 0 (with res = 0) unless 1 + x is strictly
   positive; asin and acos unless x lies in [-1, 1] (an inexact ball
   strictly inside), acosh unless x >= 1 (strictly above for an inexact
   ball) and atanh unless |x| < 1; else 1.  expm1, sinh and cosh throw
   like exp for x >= 2^(FLINT_BITS - 5) resp. |x| >= 2^(FLINT_BITS - 5). */
int mp_real_log1p_bits(mp_real_t res, const mp_real_t x, slong prec);
void mp_real_expm1_bits(mp_real_t res, const mp_real_t x, slong prec);
void mp_real_sinh_cosh_bits(mp_real_t rs, mp_real_t rc, const mp_real_t x, slong prec);
void mp_real_tanh_bits(mp_real_t res, const mp_real_t x, slong prec);
int mp_real_asin_bits(mp_real_t res, const mp_real_t x, slong prec);
int mp_real_acos_bits(mp_real_t res, const mp_real_t x, slong prec);
void mp_real_asinh_bits(mp_real_t res, const mp_real_t x, slong prec);
int mp_real_acosh_bits(mp_real_t res, const mp_real_t x, slong prec);
int mp_real_atanh_bits(mp_real_t res, const mp_real_t x, slong prec);

/* elementary functions of floating-point mantissas **************************/

/* exp(x) for x = (-1)^sgnbit (d, n) 2^(e - n FLINT_BITS) with
   e <= MP_REAL_EXP_FP_MAX_EXP = FLINT_BITS - 3: sets Y = (y, nk + 1) (y
   needs nk + 2 limbs; nk = _mp_real_elem_limbs(pb) throughout) and *err <= 9r + 104 (r the kernel's bitwise
   parameter, below 2^11 in all), and returns k such that
   exp(x) = Y 2^(k - FLINT_BITS nk) within *err units.  Y lies within
   *err units of [B^nk, 2 B^nk) (so that y[nk] <= 2) except on the paths
   with k = 0 for |x| < 1, where Y < e B^nk.  The reduction |x| = q log 2 + t uses log 2 at nk + 1
   limbs from the shared tables (exp_fp.c). */
#define MP_REAL_EXP_FP_MAX_EXP (FLINT_BITS - 3)
slong _mp_real_exp_fp(nn_ptr y, ulong * err, nn_srcptr d, slong n, slong e, int sgnbit, slong pb);

/* exp(x) for FLINT_BITS - 3 < e: exp(x) = Y 2^(K - FLINT_BITS nk) with
   the bignum K = (-1)^(returned sign) (q, *qn), q needing
   (e + FLINT_BITS - 1) / FLINT_BITS + 2 limbs; Y and *err as above (Y
   within *err units of [B^nk, 2 B^nk)).  Costs two products of about e + FLINT_BITS nk bits
   (the quotient from the tabulated 1/log 2). */
int _mp_real_exp_fp_huge(nn_ptr y, ulong * err, nn_ptr q, slong * qn, nn_srcptr d, slong n, slong e, int sgnbit, slong pb);

/* expm1(x) = (-1)^sgnbit Y 2^(k - FLINT_BITS nk) with
   Y - *elo <= |expm1(x)| 2^(FLINT_BITS nk - k) <= Y + *ehi, for
   e <= MP_REAL_EXP_FP_MAX_EXP; k = 0 for |x| < 1, where the result keeps
   its relative accuracy when the caller adds -e bits to pb */
slong _mp_real_expm1_fp(nn_ptr y, ulong * elo, ulong * ehi, nn_srcptr d, slong n, slong e, int sgnbit, slong pb);

/* 2^x = Y 2^(k - FLINT_BITS nk) within *err units, for e <= FLINT_BITS - 3;
   an integer x gives Y = B^nk, *err = 0 and k = x exactly */
slong _mp_real_exp2_fp(nn_ptr y, ulong * err, nn_srcptr d, slong n, slong e, int sgnbit, slong pb);

/* log(x) resp. log2(x) (base2) for x = (d, n) 2^(e - n FLINT_BITS) > 0,
   |e| <= 2^(FLINT_BITS - 3): sets Y = (y, *len) (y needs nk + n + 5
   limbs), *ue the radix-2 exponent of its unit and *err, and returns 1
   when log x < 0, 0 otherwise, so that
   log x = (-1)^sign Y 2^(*ue) within *err units.  x = 1 gives *len = 0.
   Near 1 the kernel runs at as many more limbs as the cancellation
   costs (at most n more: the leading zero bits of |x - 1|, which the
   callers of nearly-exact precision handle directly), so the result
   keeps its relative accuracy.  err <= 3r + 70 (r the bitwise
   parameter), times 3/2 plus 3 for log2: below 2^12 in all (log_fp.c);
   y has _mp_real_log_fp_ylen(pb, n) limbs. */
int _mp_real_log_fp(nn_ptr y, slong * len, slong * ue, ulong * err, nn_srcptr d, slong n, slong e, int base2, slong pb);

/* log(x) for x = (d, n) 2^(e - n FLINT_BITS) with a bignum e - 1 =
   (-1)^cneg (c, cn) (c normalised, |c| >= 2), as above with y needing
   _mp_real_log_fp_huge_ylen(pb, cn) = nk + cn + 3 limbs; the result has
   cn + nk limbs */
int _mp_real_log_fp_huge(nn_ptr y, slong * len, slong * ue, ulong * err, nn_srcptr d, slong n, nn_srcptr c, slong cn, int cneg, slong pb);

/* log1p(x) for x = (-1)^sgnbit (d, n) 2^(e - n FLINT_BITS) > -1 with
   e < 2^(FLINT_BITS - 3), as _mp_real_log_fp: |x| < 1/2 by the kernel of
   log(1 +- s) at -e more bits, else on 1 + x formed exactly when e <=
   FLINT_BITS nk + 2, beyond as log x plus a term below a unit; y needs
   _mp_real_log1p_fp_ylen(pb, n, e) limbs */
int _mp_real_log1p_fp(nn_ptr y, slong * len, slong * ue, ulong * err, nn_srcptr d, slong n, slong e, int sgnbit, slong pb);

/* The magnitude limit of the trigonometric functions of x (not those of
   pi x): the reduction modulo pi/2 of |x| < 2^e needs pi to about e
   bits, so mp_real_sin_cos_bits, mp_real_tan_bits and the callers of the
   _fp functions below (arf, arb, nfloat) give up (a [0 +- 1] ball, nan
   or GR_UNABLE) beyond 2^MP_REAL_TRIG_MAGLIM(prec), which keeps the cost
   polynomial in prec */
#define MP_REAL_TRIG_MAGLIM(prec) FLINT_MAX(65536, 4 * (prec))

/* Trigonometric functions of the magnitude |x| = (d, n) 2^(e - n FLINT_BITS)
   (trig_fp.c), d normalised, -WORD_MAX / 2 <= e and e <=
   MP_REAL_TRIG_MAGLIM(pb) for sin_cos and tan (the reduction costs
   O(M(e)) time; the functions of pi x reduce exactly at any e): the caller applies the sign of x
   to the odd functions.  The outputs are fixed-point values (y, *w + 1)
   at *w fraction limbs within *err ulps, nk <= *w <= wmax, where the
   kernels add limbs beyond nk for the relative accuracy of a sine
   (tangent) near its zeros as far as wmax allows (beyond, the result is
   as accurate as wmax limbs make it and the caller sees the error
   count); for sin_cos, rel says which outputs need that relative
   accuracy (bit 0 the sine, bit 1 the cosine: an output wanted only
   for a bound costs no extra limbs).  Buffers need wmax + 2 limbs (the sine and tangent of the pi
   functions max(n, nk) + 3 as well).  The return value carries the signs: bit 0 the sine
   (tangent) negative, bit 1 the cosine negative, or MP_REAL_TRIG_EXACT
   with the quadrant a in bits 0-1 when the value is exact (sine, cosine
   = 0, 1 / 1, 0 / 0, -1 / -1, 0 for a = 0, 1, 2, 3; the tangent 0 for a
   even, a pole for a odd; with MP_REAL_TRIG_ONE the tangent +-1, bit 0
   its sign), or MP_REAL_TRIG_TINY for tiny x of the pi functions, where
   the sine resp. tangent is pi x (1 -+ delta) as the product
   (y, max(n, nk) + 3) 2^(*k) with *err = 6: the sine lies in
   (Y - 1, Y + 6) units and the tangent in [Y, Y + 6); the cosine is
   B^nk at nk limbs, above the value by less than a unit.  err <= 6r + 134 for the sine and
   cosine (r the bitwise parameter), 19 (6r + 128) / 4 + 8 for the
   tangent (the quotient of the sine and cosine kernels, at most 22
   limbs on 64-bit) resp. 8r + 262 (the tangent kernel).  *recip = 1 means that tan |x| is (-1)^sign / Y instead of
   (-1)^sign Y: the kernel chooses the limbs so that Y normalises with a
   left shift of at most MP_REAL_TRIG_FP_NORM_SHIFT bits (when wmax
   allows). */
#define MP_REAL_TRIG_EXACT 8
#define MP_REAL_TRIG_TINY 16
#define MP_REAL_TRIG_ONE 32
#define MP_REAL_TRIG_FP_NORM_SHIFT (FLINT_BITS - 24)
int _mp_real_sin_cos_fp(nn_ptr ys, nn_ptr yc, slong * w, ulong * err, nn_srcptr d, slong n, slong e, slong pb, slong wmax, int rel);
int _mp_real_tan_fp(nn_ptr y, slong * w, ulong * err, int * recip, nn_srcptr d, slong n, slong e, slong pb, slong wmax);
int _mp_real_sin_cos_pi_fp(nn_ptr ys, nn_ptr yc, slong * w, slong * k, ulong * err, nn_srcptr d, slong n, slong e, slong pb, slong wmax, int rel);
int _mp_real_tan_pi_fp(nn_ptr y, slong * w, slong * k, ulong * err, int * recip, nn_srcptr d, slong n, slong e, slong pb, slong wmax);

/* Arctangents (atan_fp.c).  _mp_real_atan_fp sets (y, w + 1) to
   atan |x| B^w within err ulps for |x| = (d, n) 2^(e - n FLINT_BITS),
   with w = min(nk + limbs(-e), wmax) for |x| < 1 (the relative accuracy
   of small arguments; the caller handles tiny x) and w = nk otherwise;
   y needs wmax + 1 limbs.  _mp_real_atan2_fp evaluates theta = c pi/4 +
   sigma atan(t) for t = (Mnum / Mden) 2^-d < 1 with the mantissas
   Mnum = (num, nn), Mden = (den, dn) in [1/2, 1), d >= 0, c in {0, 2, 4}
   and sigma = +-1 (the caller sorts |x|, |y| and applies the sign of y):
   (y, len) 2^ue within elo units below and ehi above.  For c = 0 and a
   tiny t (2(d - 1) >= FLINT_BITS (nk + 1) + 4) this is the ratio itself
   at nk + 1 limbs (theta = t (1 - delta), delta t below a unit); else
   atan at w = min(nk + limbs(d - 1), wmax) (c = 0) resp. nk limbs, with
   len = w + 1, ue = -FLINT_BITS w and elo = ehi.  y needs wmax + 2 limbs
   (wmax >= nk + 1). */
void _mp_real_atan_fp(nn_ptr y, slong * w, ulong * err, nn_srcptr d, slong n, slong e, slong pb, slong wmax);
void _mp_real_atan2_fp(nn_ptr y, slong * len, slong * ue, ulong * elo, ulong * ehi, nn_srcptr num, slong nn, nn_srcptr den, slong dn, slong d, int c, int sigma, slong pb, slong wmax);

/* For (y, len) > 0 known within err units of 2^ue, sets (t, wn) to its
   leading limbs normalised (top bit set), so that Y = T 2^(*te) within
   the returned number of units of 2^(*te) (shifted-out bits count as one
   more unit).  Returns UWORD_MAX if the scaled bound would exceed
   2^(FLINT_BITS - 3). */
ulong _mp_real_fp_normalize(nn_ptr t, slong * te, nn_srcptr y, slong len, slong ue, ulong err, slong wn);

/* Hyperbolic functions of |x| = (d, n) 2^(e - n FLINT_BITS), e <=
   MP_REAL_EXP_FP_MAX_EXP, not tiny (hyperbolic_fp.c): from E = exp |x|
   and one approximate reciprocal, at w = nk limbs for |x| >= 1 and
   w = nk + ceil((1 - e) / FLINT_BITS) for |x| < 1 (the cancellation in
   E - 1/E).  _mp_real_sinh_cosh_fp sets Ys = (ys, w + 1) and
   Yc = (yc, w + 1) (either NULL; w + 2 limbs each) with sinh |x| =
   Ys 2^(k - FLINT_BITS w) and cosh x = Yc 2^(k - FLINT_BITS w) within
   *err units, returning k (-1 for |x| < 1); _mp_real_tanh_fp sets
   tanh |x| = (q, w + 1) B^-w within *err units (q needs w + 2 limbs).
   err <= 2 (that of the exponential + 3) + 2. */
slong _mp_real_sinh_cosh_fp(nn_ptr ys, nn_ptr yc, slong * w, ulong * err, nn_srcptr d, slong n, slong e, slong pb);
void _mp_real_tanh_fp(nn_ptr q, slong * w, ulong * err, nn_srcptr d, slong n, slong e, slong pb);

/* |x|^y = exp(t), t = y log |x|, for the mantissas x = (xd, xn) 2^(xe -
   xn FLINT_BITS) (|xe| <= 2^(FLINT_BITS - 3)) and y = (-1)^ysgn (yd, yn)
   2^(ye - yn FLINT_BITS) (pow_fp.c): sets Y = (y, nk + 1) (y needs nk + 2
   limbs) and *err such that |x|^y = Y 2^(K - FLINT_BITS nk) within *err
   units, with the exponent K = (-1)^(returned sign) (q, *qn), q needing
   (tbmax + FLINT_BITS - 1) / FLINT_BITS + 2 limbs (*qn = 0 for K = 0).
   Instead of the sign it returns MP_REAL_POW_ONE for |x| = 1 (the
   result is 1), MP_REAL_POW_TINY when |t| < 2^-(FLINT_BITS nk + 4), so
   that |x|^y = 1 + delta with 0 < |delta| < 2^-(FLINT_BITS nk + 3) and
   delta of the sign *tsgn (1 for negative), or MP_REAL_POW_HUGE when
   |t| >= 2^tbmax (tbmax >= FLINT_BITS - 3) up to the accuracy of the
   computed t (a true |t| slightly below 2^tbmax may give HUGE; a cheap
   lower bound on |t| from the exponents and the leading bits of x
   decides most such arguments without evaluating log |x|), the sign of
   t in *tsgn.  The logarithm runs at nk + ceil((tb + 12) /
   FLINT_BITS) limbs for the bound 2^tb on |t|; *err is that of the
   exponential plus the propagated error of t, which can reach about
   2^57 units when the units are finer than pb bits (the relative error
   stays below 2^-(pb - 12)). */
#define MP_REAL_POW_ONE 2
#define MP_REAL_POW_TINY 3
#define MP_REAL_POW_HUGE 4
int _mp_real_pow_fp(nn_ptr y, ulong * err, nn_ptr q, slong * qn, int * tsgn, nn_srcptr xd, slong xn, slong xe, nn_srcptr yd, slong yn, slong ye, int ysgn, slong pb, slong tbmax);

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

/* sin(pi p/q) and cos(pi p/q) for 0 <= p <= q/2: outputs (ys, n + 1),
   (yc, n + 1), either may be NULL.  err <= 6r + 130, r <= 768 the
   bitwise parameter, when the sin/cos kernel is used, else at most 4
   (computed at run time).  alg: 0 = tuned, 1 = the kernel on pi p/q,
   2 = the Chebyshev (U_(q-1) root) iteration of order r (0 =
   default), 3 = square roots for q = 5, 8, 10, 12. */
void _mp_real_sin_cos_pi_ui_div_ui(nn_ptr ys, nn_ptr yc, ulong * err, ulong p, ulong q, slong n);
void _mp_real_sin_cos_pi_ui_div_ui_tune(nn_ptr ys, nn_ptr yc, ulong * err, ulong p, ulong q, slong n, int alg, int r);

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
   halving depth r >= 16 (tuning). */
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

/* the logarithms log p_j of the first num primes and the angles pi/2,
   2 arg(pi_j) of the first num nonreal Gaussian primes, j < num, as arb
   balls at prec bits (arb_log_primes_vec_bsplit and
   arb_atan_gauss_primes_vec_bsplit wrap these) */
void _mp_real_log_primes_vec_arb(arb_ptr res, slong num, slong prec);
void _mp_real_atan_gauss_vec_arb(arb_ptr res, slong num, slong prec);

#ifdef __cplusplus
}
#endif

#endif
