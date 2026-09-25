/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef FIXED_H
#define FIXED_H

#include <stdint.h>
#include "flint.h"
#include "arb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Efficient low-level fixed-point real arithmetic.

   A fixed-point number (x, n) is an unsigned n-limb fraction
   x[0], ..., x[n-1] representing sum x[i] 2^(FLINT_BITS (i - n)),
   i.e. 0 <= x < 1 with unit in the last place (ulp)
   2^(-FLINT_BITS n).  Outputs of size n + 1 additionally carry an
   integer (units) limb at index n.

   The series evaluation functions below require an argument reduced
   below 2^-32 (checked with a FLINT_ASSERT) and dispatch internally
   on the top limb of x: nonzero selects the 32-bit reduction range,
   zero the 64-bit-or-higher one.  Error bounds are in ulp; for exp
   and the non-alternating (hyperbolic) functions they are one-sided
   (the result never exceeds the true value), for the alternating
   functions two-sided. */

#define FIXED_EXP_RS_MAX_ERR(n) 10
/* for r = 0 (tuned default) the bound holds with the selected r,
   which never exceeds 512 */
#define FIXED_EXP_BITWISE_RS_MAX_ERR(n, r) \
    (9 * (slong) ((r) ? (r) : 512) + 100)
#define FIXED_LOG1P_BITWISE_RS_MAX_ERR(n, r) \
    (3 * (slong) ((r) ? (r) : 512) + 64)
#define FIXED_SIN_COS_BITWISE_RS_MAX_ERR(n, r) \
    (6 * (slong) ((r) ? (r) : 512) + 128)
#define FIXED_ATAN_BITWISE_RS_MAX_ERR(n, r) \
    (4 * (slong) ((r) ? (r) : 512) + 64)
#define FIXED_TAN_BITWISE_RS_MAX_ERR(n, r) \
    (8 * (slong) ((r) ? (r) : 512) + 256)
#define FIXED_SIN_RS_MAX_ERR(n) 15
#define FIXED_COS_RS_MAX_ERR(n) 15
#define FIXED_SIN_COS_RS_MAX_ERR(n) 15
#define FIXED_SINH_RS_MAX_ERR(n) 15
#define FIXED_COSH_RS_MAX_ERR(n) 15
#define FIXED_SINH_COSH_RS_MAX_ERR(n) 15
#define FIXED_ATAN_RS_MAX_ERR(n) 15
#define FIXED_ATANH_RS_MAX_ERR(n) 15

/* exp((x, n)) -> (res, n + 1) */
void fixed_exp_rs(nn_ptr res, nn_srcptr x, slong n);

/* exp((x, n)) -> (res, n + 1) for any 0 <= x < 1 (n >= 2 on 32-bit
   limbs; r = 0 selects a tuned default), using bitwise
   argument reduction with a runtime-cached table of log(1 + 2^-i)
   followed by Taylor evaluation below 2^-r (r >= 32 a tuning
   parameter) and shift-and-add reconstruction.  All work happens at
   the output precision, so the error bound grows linearly with r;
   callers wanting sub-ulp accuracy should pad the precision by one
   limb themselves. */
/* exp(t) of a reduced argument t < 2^-r, r >= 32, into
   (y, wn + 1): wn fraction limbs and a units limb.  alg: 0 = tuned
   automatic choice, 1 = direct rectangular-splitting series,
   2 = sinh series + square root, 3 = one bit-burst step + sinh,
   4 = full bit-burst; the arithmetic around the series kernels and
   the exact binary splitting is fball ball arithmetic, with a
   rigorous bound checked against FIXED_EXP_REDUCED_MAX_ERR. */
#define FIXED_EXP_REDUCED_MAX_ERR 96
void fixed_exp_reduced(nn_ptr y, nn_srcptr t, slong wn,
    flint_bitcnt_t r, int alg);

/* sin(t) and g = 1 - cos(t) of a reduced argument t < 2^-r,
   r >= 16 (the series algorithms 1 and 2 require r >= 32), into
   (ysin, wn) and (yg, wn); the trigonometric analogue of
   fixed_exp_reduced.  alg: 0 = tuned automatic choice, 1 = direct
   sine + cosine rectangular-splitting series, 2 = sine series plus
   a squaring and a square root, 3 = one bit-burst step + series,
   4 = full bit-burst, in fball arithmetic as fixed_exp_reduced.
   Error at most FIXED_SIN_COS_REDUCED_MAX_ERR ulps on each output. */
#define FIXED_SIN_COS_REDUCED_MAX_ERR 96
void fixed_sin_cos_reduced(nn_ptr ysin, nn_ptr yg, nn_srcptr t,
    slong wn, flint_bitcnt_t r, int alg);

/* exp, sine and cosine on [0, 1) without table-based argument
   reduction: leading-zero inspection, repeated halving to a tuned
   depth r(n), fixed_*_reduced, and squarings back (for the
   trigonometric pair, cosine doublings and one final square root
   for the sine), with internal guard limbs keeping the doubling
   amplification below one output ulp.  Outputs carry n fraction
   limbs and a unit limb.  Errors at most the *_MAX_ERR ulps. */
#define FIXED_EXP_NOTAB_MAX_ERR 128
void fixed_exp_notab(nn_ptr y, nn_srcptr x, slong n);
#define FIXED_SIN_COS_NOTAB_MAX_ERR 128
void fixed_sin_cos_notab(nn_ptr ysin, nn_ptr ycos, nn_srcptr x,
    slong n);

/* exp on [0, 1) by diophantine (multi-prime) argument reduction (the fixed-point
   port of arb_exp_arf_log_reduction, see exp_diophantine.c):
   exp(x) = 2^c_0 (p/q) exp(t) with p, q products of the first
   num_primes primes and t reduced to a tiny residual through a
   table of integer relations between their logarithms.  Output: n
   fraction limbs and a unit limb, within
   FIXED_EXP_DIOPHANTINE_MAX_ERR ulps.  The tunable worker takes the
   number of primes (2 <= num_primes <= FIXED_LOG_PRIMES_MAX; the
   relation table comes from fixed_rel_table, precomputed for the
   common counts) and the weight budget of the reduction
   (the bound on sum_{j > 0} |c_j| log2(p_j) / log 2, a proxy for
   the bit size of p q); the public entry uses 13 primes and
   FLINT_BITS n, as arb does.  Precomputation is one logarithm per
   prime at the working precision, cached per thread. */
#define FIXED_EXP_DIOPHANTINE_MAX_ERR 3
void fixed_exp_diophantine(nn_ptr y, nn_srcptr x, slong n);
void _fixed_exp_diophantine_tune(nn_ptr y, nn_srcptr x, slong n,
    slong num_primes, double max_weight);

/* Relation tables for the diophantine reductions (rel_tab.c): for
   the first num primes (gaussian = 0: alpha_j = log p_j) or the
   first num nonreal Gaussian primes (gaussian = 1: alpha_0 = pi/2,
   alpha_j = 2 arg pi_j, with pi_j = a + b i from
   _fixed_gaussian_primes), the rows of d (num entries each) are
   integer relations sum_j d_ij alpha_j = epsilon_i with |epsilon_i|
   decreasing.  Tables for num = 2, 4, 6, 8, 10, 12, 13, 16, 20, 24,
   32, 40, 48 are precomputed; other sizes are generated on first
   use (seconds around 48 primes) and cached per thread, which
   fixed_rel_table_is_cached predicts. */
#define FIXED_REL_MAX 64
#define FIXED_LOG_PRIMES_MAX FIXED_REL_MAX
#define FIXED_REL_TERMINATOR -32768
typedef struct
{
    slong num;
    slong rows;
    int gaussian;
    int is_static;
    const ulong * primes;        /* the rational primes; NULL if gaussian */
    const float * weights;       /* log2 p_j, or log N(pi_j); weights[0] = 0 */
    const short * d;             /* rows x num, then a row starting with
                                    FIXED_REL_TERMINATOR */
    const double * epsilon;
    const double * epsilon_inv;
    double epsilon_min;          /* |epsilon| of the last row */
}
fixed_rel_struct;

const fixed_rel_struct * fixed_rel_table(int gaussian, slong num);
int fixed_rel_table_is_cached(int gaussian, slong num);

/* Machin-type sets (machin_tab.c): log p_i = (1/den) sum_j c[i][j]
   atanh(1/x_j) for the first num primes (gaussian = 0), or
   arg pi_i = (1/den) sum_j c[i][j] atan(1/x_j) for the first num
   nonreal Gaussian primes (gaussian = 1); c is num x num.
   fixed_machin_table returns the best set for num values: the largest
   one with at most num terms (the smallest set when num is below
   that), the remaining values being left to one followup series each.
   Sets exist for every num from 4 (3 for the Gaussian primes) to 32,
   and for 40 and 48; fixed_machin_table_max gives the largest. */
#define FIXED_MACHIN_X_LIMBS (128 / FLINT_BITS)
#define FIXED_MACHIN_MAX_NUM 64
#define FIXED_MACHIN_MAX_PRIMES 8
typedef struct
{
    slong num;
    const ulong * x;             /* the arguments, FIXED_MACHIN_X_LIMBS
                                    limbs (128 bits) each */
    ulong den;
    int cbits;                   /* one more than the largest bit length
                                    among the coefficients, which are
                                    reconstructed on first use */
    int gaussian;
}
fixed_machin_struct;

const fixed_machin_struct * fixed_machin_table(int gaussian, slong num);
slong fixed_machin_table_max(int gaussian);
/* the arguments and coefficients as fmpz, whatever their size; the
   row function fills row[0..num) with c[i][0..num) */
void fixed_machin_get_x(fmpz_t q, const fixed_machin_struct * tab, slong i);
void fixed_machin_get_c(fmpz_t c, const fixed_machin_struct * tab, slong i,
    slong j);
void fixed_machin_get_c_row(fmpz * row, const fixed_machin_struct * tab,
    slong i);
nn_srcptr fixed_machin_c_row_raw(const fixed_machin_struct * tab, slong i,
    slong * climbs, const unsigned char ** csign);

/* the precomputed tables (rel_tab_data.c) */
typedef struct
{
    int gaussian;
    slong num;
    slong rows;
    const short * d;
    const double * epsilon;
}
fixed_rel_static_struct;
FLINT_DLL extern const fixed_rel_static_struct _fixed_rel_static[];
FLINT_DLL extern const slong _fixed_rel_static_num;

/* real and imaginary parts, consecutively, of the first 64 nonreal
   Gaussian primes in order of norm: 1+i, 1+2i, 2+3i, ... */
#define FIXED_ATAN_GAUSS_MAX 64
/* FLINT_DLL as for the other exported tables: the build exports
   functions automatically on Windows, but not data symbols */
FLINT_DLL extern const signed char
    _fixed_gaussian_primes[2 * FIXED_ATAN_GAUSS_MAX];

/* Internal: thread-local cache of the angles pi/2, 2 arg(pi_j)
   (atan_gauss.c), laid out like the logarithm cache below;
   _fixed_atan_gauss_vec gives the angles as arb balls (from
   _fixed_atan_gauss_vec_fball: the Machin-type sets of machin_tab.c
   for up to 48, single arctangents beyond). */
void _fixed_atan_gauss_vec(arb_ptr res, slong num, slong prec);
void _fixed_atan_gauss_ensure(slong num, slong nv);
nn_srcptr _fixed_atan_gauss_entry(slong j, slong nv);
void _fixed_atan_gauss_clear(void);

/* sin and cos of (x, n) in [0, 1) by diophantine (multi-prime)
   argument reduction (sin_cos_diophantine.c): x = c_0 pi/2 +
   sum_j c_j 2 arg(pi_j) + t with t tiny, and e^(ix) = i^c_0 e^(it)
   A^2 / |A|^2 for the Gaussian integer A = prod pi_j^(c_j)
   (conjugates for c_j < 0), whose norm is a rational integer.
   Outputs (either may be NULL) carry n fraction limbs and a unit
   limb, within FIXED_SIN_COS_DIOPHANTINE_MAX_ERR ulps.  The tunable
   workers take the number of Gaussian primes and the weight budget
   (arb uses 13 and half the precision; the default here is 32 and
   four times the precision). */
#define FIXED_SIN_COS_DIOPHANTINE_MAX_ERR 4
void fixed_sin_cos_diophantine(nn_ptr ysin, nn_ptr ycos, nn_srcptr x,
    slong n);
void _fixed_sin_cos_diophantine_tune(nn_ptr ysin, nn_ptr ycos,
    nn_srcptr x, slong n, slong num_primes, double max_weight);

/* tan of (x, n) in [0, 1) -> (res, n + 1) by the same reduction: the
   ratio of the two parts of e^(it) A^2, one division and no
   normalization at all */
void fixed_tan_diophantine(nn_ptr res, nn_srcptr x, slong n);
void _fixed_tan_diophantine_tune(nn_ptr res, nn_srcptr x, slong n,
    slong num_primes, double max_weight);

/* Newton-Taylor inverses from the forward functions (newton.c):
   -log(x) for (x, n) in [1/2, 1) -> (y, n) and atan(x) for (x, n) in
   [0, 1) -> (y, n), each one step of order N from a starting value
   at about 1/(N + 1) (log) resp. 1/(2N) (atan) of the precision --
   the bitwise function up to FIXED_NEWTON_CUTOFF limbs, the step
   itself recursively above -- followed by one forward evaluation
   (exp resp. sin_cos) at the working precision of n + 1 limbs and
   a short Taylor polynomial in the residual, all in fball
   arithmetic with a rigorous bound checked on export.  Errors at
   most the *_MAX_ERR ulps.  The tunable workers take the forward
   algorithm (0 = diophantine, 1 = bitwise, 2 = notab) and the number N of
   series terms (0 = default). */
#define FIXED_NEGLOG_NEWTON_MAX_ERR 2
#define FIXED_ATAN_NEWTON_MAX_ERR 2
#ifndef FIXED_NEWTON_CUTOFF
#define FIXED_NEWTON_CUTOFF 512
#endif
void fixed_neglog_newton(nn_ptr y, nn_srcptr x, slong n);
/* -log(x) for (x, n) in [1/2, 1) by the Sasaki-Kanada formula
   (log_agm.c): one AGM of theta functions, pi/4 and log 2 from the
   per-thread caches; within FIXED_NEGLOG_NEWTON_MAX_ERR ulps.  The
   tunable worker takes the number N of theta_3 terms (0 = default) */
void fixed_neglog_agm(nn_ptr y, nn_srcptr x, slong n);
void _fixed_neglog_agm_tune(nn_ptr y, nn_srcptr x, slong n, slong N);
void fixed_atan_newton(nn_ptr y, nn_srcptr x, slong n);
void _fixed_neglog_newton_tune(nn_ptr y, nn_srcptr x, slong n, int forward,
    slong N);
void _fixed_atan_newton_tune(nn_ptr y, nn_srcptr x, slong n, int forward,
    slong N);

/* Internal: the relation-table descent shared by the diophantine (multi-prime)
   reductions (exp_diophantine.c); the angles (log p_j, or 2 arg of
   the Gaussian primes) are read as wr + 1 limb entries at
   alpha + j stride, unit limb on top. */
slong _fixed_log_reduce(slong * rel, const fixed_rel_struct * tab,
    nn_srcptr x, slong wr, double max_weight, double eps_min,
    nn_srcptr alpha, slong stride);
void _fixed_log_dot(nn_ptr acc, nn_srcptr base, slong len,
    const slong * rel, slong num, nn_srcptr alpha, slong stride);
double _fixed_signed_get_d(nn_srcptr a, slong len, nn_ptr tmp);

/* Internal: thread-local cache of fixed-point logarithms of the
   first primes 2, 3, 5, ... (log_primes.c).  _ensure(num, nv) makes
   the table cover num primes with at least nv fraction value limbs
   each (plus a guard limb below and a unit limb above);
   _entry(j, nv) returns the top nv + 1 limbs of entry j -- nv
   fraction limbs with the unit limb at index nv -- exactly
   floor(log(p_j) B^nv) spread over them, valid until the next
   _ensure call on this thread. */
void _fixed_log_primes_ensure(slong num, slong nv);

nn_srcptr _fixed_log_primes_entry(slong j, slong nv);
void _fixed_log_primes_clear(void);
slong _fixed_log_primes_max_limbs(void);

/* internal forced-depth workers (the public entry points choose r
   from tuned tables; these take it explicitly, for tuning) */
void _fixed_exp_notab_r(nn_ptr y, nn_srcptr x, slong n, int r);
void _fixed_sin_cos_notab_r(nn_ptr ysin, nn_ptr ycos, nn_srcptr x,
    slong n, int r);

/* notab tuning: the trigonometric bit-burst boundary and the sine
   square root's Newton cutoff */
#ifndef FIXED_SIN_COS_NOTAB_BURST_CUTOFF
#define FIXED_SIN_COS_NOTAB_BURST_CUTOFF 1200
#endif
#ifndef FIXED_SIN_COS_NOTAB_SQRT_NEWTON_CUTOFF
#define FIXED_SIN_COS_NOTAB_SQRT_NEWTON_CUTOFF 2000
#endif

void fixed_exp_bitwise_rs(nn_ptr res, nn_srcptr x, slong n, int r);

/* fully specialized per-size implementations (default dispatch;
   generated and tuned by dev/tune_fixed.py) */
#if FLINT_BITS == 64
void fixed_exp_opt_1(nn_ptr res, nn_srcptr x);
void fixed_exp_opt_2(nn_ptr res, nn_srcptr x);
void fixed_exp_opt_3(nn_ptr res, nn_srcptr x);
void fixed_exp_opt_4(nn_ptr res, nn_srcptr x);
void fixed_exp_opt_5(nn_ptr res, nn_srcptr x);
void fixed_exp_opt_6(nn_ptr res, nn_srcptr x);
void fixed_exp_opt_7(nn_ptr res, nn_srcptr x);
void fixed_atan_opt_1(nn_ptr res, nn_srcptr x);
void fixed_atan_opt_2(nn_ptr res, nn_srcptr x);
void fixed_atan_opt_3(nn_ptr res, nn_srcptr x);
void fixed_atan_opt_4(nn_ptr res, nn_srcptr x);
void fixed_atan_opt_5(nn_ptr res, nn_srcptr x);
void fixed_atan_opt_6(nn_ptr res, nn_srcptr x);
void fixed_atan_opt_7(nn_ptr res, nn_srcptr x);
void fixed_log1p_opt_1(nn_ptr res, nn_srcptr x);
void fixed_log1p_opt_2(nn_ptr res, nn_srcptr x);
void fixed_log1p_opt_3(nn_ptr res, nn_srcptr x);
void fixed_log1p_opt_4(nn_ptr res, nn_srcptr x);
void fixed_log1p_opt_5(nn_ptr res, nn_srcptr x);
void fixed_log1p_opt_6(nn_ptr res, nn_srcptr x);
void fixed_log1p_opt_7(nn_ptr res, nn_srcptr x);
void _fixed_trig_opt_1(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x);
void _fixed_trig_opt_2(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x);
void _fixed_trig_opt_3(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x);
void _fixed_trig_opt_4(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x);
void _fixed_trig_opt_5(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x);
void _fixed_trig_opt_6(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x);
void _fixed_trig_opt_7(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x);
void _fixed_trig_opt_8(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x);
void _fixed_trig_opt_9(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x);
void _fixed_trig_opt_10(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x);
void _fixed_trig_opt_11(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x);
void _fixed_trig_opt_12(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x);
void fixed_sin_cos_opt_1(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void fixed_sin_cos_opt_2(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void fixed_sin_cos_opt_3(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void fixed_sin_cos_opt_4(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void fixed_sin_cos_opt_5(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void fixed_sin_cos_opt_6(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void fixed_sin_cos_opt_7(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void fixed_sin_cos_opt_8(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void fixed_sin_cos_opt_9(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void fixed_sin_cos_opt_10(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void fixed_sin_cos_opt_11(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void fixed_sin_cos_opt_12(nn_ptr ysin, nn_ptr ycos, nn_srcptr x);
void fixed_tan_opt_1(nn_ptr res, nn_srcptr x);
void fixed_tan_opt_2(nn_ptr res, nn_srcptr x);
void fixed_tan_opt_3(nn_ptr res, nn_srcptr x);
void fixed_tan_opt_4(nn_ptr res, nn_srcptr x);
void fixed_tan_opt_5(nn_ptr res, nn_srcptr x);
void fixed_tan_opt_6(nn_ptr res, nn_srcptr x);
void fixed_tan_opt_7(nn_ptr res, nn_srcptr x);
void fixed_tan_opt_8(nn_ptr res, nn_srcptr x);
void fixed_tan_opt_9(nn_ptr res, nn_srcptr x);
void fixed_tan_opt_10(nn_ptr res, nn_srcptr x);
void fixed_tan_opt_11(nn_ptr res, nn_srcptr x);
void fixed_tan_opt_12(nn_ptr res, nn_srcptr x);
#endif

/* shared internals of the specialized per-size implementations */
void _fixed_exp_recon(nn_ptr y, nn_ptr sh, slong ylen, const slong * used,
    slong j, slong num);
void _fixed_tan_halfangle_mid(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x, slong n, int r, void (*series)(nn_ptr, nn_srcptr));

/* the reduction parameter that r = 0 selects at size n: the
   compile-time constant of the specialized per-size implementation
   where one exists, the tuned large-n ladder beyond (multiples of 64
   chosen by src/fixed/tune/tune-bitwise-r.c) */
int fixed_exp_bitwise_rs_default_r(slong n);
int fixed_log1p_bitwise_rs_default_r(slong n);
int fixed_atan_bitwise_rs_default_r(slong n);
int fixed_trig_bitwise_rs_default_r(slong n);

/* log(1 + (x, n)) -> (res, n) for any 0 <= x < 1, by the dual
   reduction: greedily multiply P by factors 1 + 2^-i (each one
   shift-and-add) while P (1 + 2^-i) <= 1 + x, then
   log((1+x)/P) = 2 atanh(((1+x) - P)/((1+x) + P)) with a single
   division, then add the tabulated logarithms.  The error bound
   grows linearly with r as for fixed_exp_bitwise_rs.  Requires
   r = 0 (which selects a tuned default) or r >= 16; values below 32 shorten the reduction further and are
   effective in the specialized code for n <= 4. */
void fixed_log1p_bitwise_rs(nn_ptr res, nn_srcptr x, slong n, int r);

/* sin and cos of (x, n) in [0, 1) -> (ysin, n + 1), (ycos, n + 1)
   (either may be NULL), by greedy reduction with the angles
   atan(2^-i) and the tangent half-angle reconstruction (see
   tan_bitwise_rs.c); requires n >= 2 on 32-bit limbs and r = 0
   (tuned default) or r >= 16. */
void fixed_sin_cos_bitwise_rs(nn_ptr ysin, nn_ptr ycos, nn_srcptr x,
    slong n, int r);

/* atan((x, n)) -> (res, n) for x in [0, 1), by greedy vectoring;
   r = 0 selects a tuned default; otherwise r >= 16. */
void fixed_atan_bitwise_rs(nn_ptr res, nn_srcptr x, slong n, int r);

/* tan((x, n)) -> (res, n + 1) for x in [0, 1); tan(1) < 1.56, so the
   result carries a unit limb.  r = 0 selects a tuned default. */
void fixed_tan_bitwise_rs(nn_ptr res, nn_srcptr x, slong n, int r);

/* Internal: sin, cos and tan of (x, n) in [0, 1) by the tangent
   half-angle reconstruction; any output may be NULL.  Returns 1 if the
   size is handled, 0 if the caller must fall back.  Available on both
   word sizes: off 64-bit limbs there is no hand-written tangent series,
   and tan(t') comes from the sine and cosine series instead. */
int _fixed_tan_halfangle(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan,
    nn_srcptr x, slong n, int r);

/* sin, cos, sinh, cosh of (x, n) -> (res, n + 1); the combined
   versions allow either output to be NULL */
void fixed_sin_rs(nn_ptr res, nn_srcptr x, slong n);
void fixed_cos_rs(nn_ptr res, nn_srcptr x, slong n);
void fixed_sin_cos_rs(nn_ptr ysin, nn_ptr ycos, nn_srcptr x, slong n);
void fixed_sinh_rs(nn_ptr res, nn_srcptr x, slong n);
void fixed_cosh_rs(nn_ptr res, nn_srcptr x, slong n);
void fixed_sinh_cosh_rs(nn_ptr ysinh, nn_ptr ycosh, nn_srcptr x, slong n);

/* atan, atanh of (x, n) -> (res, n) */
void fixed_atan_rs(nn_ptr res, nn_srcptr x, slong n);
void fixed_atanh_rs(nn_ptr res, nn_srcptr x, slong n);

/* Fallbacks used on all architectures and precisions: rectangular
   splitting at constant full precision with coefficients generated on
   the fly, requiring only x < 2^-32 (the number of terms is chosen
   from the actual leading zero bits of x).  Exposed for testing. */
void _fixed_exp_rs_fallback(nn_ptr res, nn_srcptr x, slong n);

/* Internal: exp for the wider range x < 2^-16 (n <= 5), used by the
   small reduction parameters of fixed_exp_bitwise_rs.  The hardcoded
   series family exists only for 64-bit limbs. */
#if FLINT_BITS == 64

/* Internal: the fully specialized exp series, one per n <= 5, each
   built for the reduction parameter hardcoded alongside it in
   fixed_exp_bitwise_rs and using the smallest number of terms N with
   N r + log2(N!) >= 64 n. */
#endif

/* Internal: thread-local cached table of L_i = log(1 + 2^-i), one
   entry per index i = 0..r, shared by fixed_exp_bitwise_rs and
   fixed_log1p_bitwise_rs.  _ensure(nv, r) makes the table cover the
   indices 0..r with at least nv value limbs each (plus a guard limb
   below them).

   The storage itself is thread-local and is deliberately NOT
   declared here: Windows DLLs cannot export thread-local data, so
   library-external code (the test suite) reads entries through the
   accessors below, while the module's own translation units see the
   definitions via src/fixed/impl.h.  _entry(i, n) returns
   the top n limbs of entry i, valid until the next _ensure call on
   this thread. */
void _fixed_exp_logs_ensure(slong nv, slong rc);
nn_srcptr _fixed_exp_logs_entry(slong i, slong n);
void _fixed_exp_logs_clear(void);
slong _fixed_exp_logs_max_index(void);

/* Internal: number of slots the used array of _fixed_bitwise_reduce
   must provide: each index i = istart..r is used at most once, plus
   the window-boundary and final steps, which may repeat an index a
   bounded number of times to absorb the truncation creep of the
   table (see exp_bitwise_rs.c). */
#define FIXED_BITWISE_REDUCE_USED_ALLOC(r) \
    ((r) + 2 * ((r) / FLINT_BITS) + 12)

/* Internal: shared greedy table-subtraction reduction (see
   exp_bitwise_rs.c); returns the number of indices recorded in
   used, which must have room for
   FIXED_BITWISE_REDUCE_USED_ALLOC(r) entries. */
slong _fixed_bitwise_reduce(nn_ptr t, slong wn, int r, slong istart,
    nn_srcptr tab, slong tabn, slong * used);

/* Internal: thread-local cached table of the angles
   A_i = atan(2^-i) (entry 0, A_0 = pi/4, is unused by the
   reductions, which start at i = 1, but fits the fraction format
   and is tabulated anyway), shared by fixed_sin_cos_bitwise_rs,
   fixed_tan_bitwise_rs and fixed_atan_bitwise_rs.  Storage and
   accessors work exactly as for the logarithm table above. */
void _fixed_atans_ensure(slong nv, slong rc);

/* Internal: n-limb one-sided fixed-point approximations of the table
   values by binary splitting in fball arithmetic (floor(v B^n) or one
   below it, never above); i >= 1. */
void fixed_atan_2mexp_ui_bs(nn_ptr res, ulong i, slong n);
void fixed_log1p_2mexp_ui_bs(nn_ptr res, ulong i, slong n);

/* Approximate (not ulp-accurate) fixed-point inversion,
   division and square roots by Newton / Karp-Markstein iteration on
   middle products; ports of the radix_*_approx functions.  Writing
   B = 2^64:

   fixed_inv_newton: given (a, an), a_{an-1} != 0, representing
   a in [1/B, 1) with an fraction limbs, sets (q, n+2) to 1/a in
   (1, B] with n fraction limbs and two integral limbs (the top limb
   may be zero); |error| <= 4 B^-n / a.

   fixed_div_newton: numerator (b, bn) in [0, 1), denominator
   (a, an) in [1/B, 1) with a_{an-1} != 0; sets (q, n+2) to b/a with
   n fraction limbs and two integral limbs; |error| <= 4 B^-n / a.

   fixed_rsqrt_ui_newton: 2 <= a < B; sets (res, n) to the fraction
   limbs of 1/sqrt(a); |error| <= 2 B^-n.

   fixed_rsqrt_newton: (a, an) in [B^-2, 1), one of the two top limbs
   nonzero; sets (q, n+2) to 1/sqrt(a) in (1, B] with n fraction
   limbs and two integral limbs; |error| <= 4 B^-n / sqrt(a).

   fixed_sqrt_newton: input as for fixed_rsqrt_newton; sets (q, n+2)
   to sqrt(a) in [1/B, 1) (the value can round up to 1);
   |error| <= 4 B^-n / sqrt(a).
*/
void fixed_inv_newton_basecase(nn_ptr q, nn_srcptr a, slong an, slong n);
void fixed_inv_newton(nn_ptr q, nn_srcptr a, slong an, slong n);
void fixed_div_newton_invmul(nn_ptr q, nn_srcptr b, slong bn, nn_srcptr a, slong an, slong n);
void fixed_div_newton(nn_ptr q, nn_srcptr b, slong bn, nn_srcptr a, slong an, slong n);
void fixed_rsqrt_ui_newton_basecase(nn_ptr res, ulong a, slong n);
void fixed_rsqrt_ui_newton(nn_ptr res, ulong a, slong n);
void fixed_rsqrt_newton_basecase(nn_ptr q, nn_srcptr a, slong an, slong n);
void fixed_rsqrt_newton(nn_ptr q, nn_srcptr a, slong an, slong n);
void fixed_sqrt_newton_rsqrtmul(nn_ptr q, nn_srcptr a, slong an, slong n);
void fixed_sqrt_newton(nn_ptr q, nn_srcptr a, slong an, slong n);
nn_srcptr _fixed_atans_entry(slong i, slong n);
void _fixed_atans_clear(void);
slong _fixed_atans_max_index(void);
void _fixed_sin_cos_rs_fallback(nn_ptr ysin, nn_ptr ycos, nn_srcptr x,
    slong n, int alternating);
void _fixed_atan_rs_fallback(nn_ptr res, nn_srcptr x, slong n,
    int alternating);

/* Internal: atanh for the wider range x < 2^-16, used by the small
   reduction parameters of fixed_log1p_bitwise_rs. */

/* Internal: atan and sin/cos for the wider range x < 2^-16, used by
   the small reduction parameters of the bitwise trigonometric
   functions.  The sin/cos routine requires both outputs. */

/* Internal: hand-written atan series, one per n <= 4, each built for
   the reduction parameter hardcoded alongside it in
   fixed_atan_bitwise_rs (64-bit limbs only). */
#if FLINT_BITS == 64

/* Internal: hand-written atanh and sin/cos series, each built for the
   reduction parameter hardcoded alongside it in the bitwise callers.
   The sin/cos routines compute both outputs from a single squaring. */

/* Internal: tan series for the half-angle reconstruction, one per
   n <= 12, each built for the reduction parameter hardcoded alongside
   it in tan_bitwise_rs.c. */
#endif


/* -------------------------------------------------------------- */
/* Internal declarations (formerly fixed/impl.h), placed here so
   that the test, tune and profile subdirectories build without
   extra include paths. */

/* Library-internal view of the cached reduction tables.

   These thread-local objects are shared across the translation units
   of the module (the dispatch files, the reduction, and the
   specialized per-size implementations) but are deliberately NOT
   declared in fixed.h: Windows DLLs cannot export thread-local data,
   so external consumers -- the test suite -- go through the
   _fixed_exp_logs_entry / _fixed_atans_entry accessors instead.

   Each entry occupies _fixed_{exp_logs,atans}_n limbs: the value
   limbs with one guard limb below them.  Consumers wanting the top n
   limbs of entry i read tab + i * stride + (stride - n). */

#define FIXED_STATIC_TAB_INLINE static inline

#ifdef __cplusplus
extern "C" {
#endif

extern FLINT_TLS_PREFIX nn_ptr _fixed_exp_logs;
extern FLINT_TLS_PREFIX slong _fixed_exp_logs_n;
extern FLINT_TLS_PREFIX slong _fixed_exp_logs_r;

extern FLINT_TLS_PREFIX nn_ptr _fixed_atans;
extern FLINT_TLS_PREFIX slong _fixed_atans_n;
extern FLINT_TLS_PREFIX slong _fixed_atans_r;

extern FLINT_TLS_PREFIX nn_ptr _fixed_log_primes;
extern FLINT_TLS_PREFIX slong _fixed_log_primes_n;
extern FLINT_TLS_PREFIX slong _fixed_log_primes_num;

extern FLINT_TLS_PREFIX nn_ptr _fixed_atan_gauss;
extern FLINT_TLS_PREFIX slong _fixed_atan_gauss_n;
extern FLINT_TLS_PREFIX slong _fixed_atan_gauss_num;

/* Static prefixes of the two tables covering all reductions with
   r <= FIXED_STATIC_TAB_R whose per-entry reads fit in
   FIXED_STATIC_TAB_N limbs.  Unlike the dynamic tables, whose
   bottom limb is a guard used for in-place generation, every stored
   limb here is a value limb (the entries are the top limbs of a
   dynamic table built one limb deeper).  The accessors below hand
   out the static data when it suffices -- avoiding the
   precomputation, the TLS lookup, and the wide entry stride of a
   high-precision dynamic table -- and fall back to ensuring the
   dynamic one. */

#define FIXED_STATIC_TAB_N 12
#define FIXED_STATIC_TAB_R 32

/* internal: mpn binary splitting for sum_{k=1}^N x^k / (k! 2^(rk))
   (exp_sum_bs.c); T needs (N (r + 128))/64 + 4 limbs, Q needs
   (N bits(N+1))/64 + 3 */
slong _fixed_exp_bs_num_terms(flint_bitcnt_t r, slong prec);
void _fixed_sin_cos_sum_bs_powtab(nn_ptr A, slong * an, slong * ae,
    nn_ptr B, slong * bn, slong * be, nn_ptr Q, slong * qn,
    slong * QE, nn_srcptr xp, slong xn, slong D, slong N,
    slong lmax);
void _fixed_exp_sum_bs_powtab(nn_ptr T, slong * tn, nn_ptr Q,
    slong * qn, slong * QE, nn_srcptr xp, slong xn, slong D,
    slong N);

/* fball: SEMI-PRIVATE ball arithmetic (subject to future revision).
   The interface below is documented for use inside FLINT and by
   expert callers, but makes no stability promises across releases:
   the representation, the error-normalization policy and the set of
   operations may all change.

   Internal mpf-like floating-point numbers with arb-like ball
   semantics, in radix B = 2^FLINT_BITS, intended as a backend for
   algorithms on scaled or growing numbers (AGM iterations, binary
   splitting summation) where all cheap operations should stay in mpn
   arithmetic and exponents are limb-aligned so that shifts are limb
   copies.

   Representation:

       value = (-1)^negative * (d, size) * B^(exp - size)

   with the mantissa (d, size) an unsigned little-endian integer whose
   top limb is nonzero, so B^(exp-1) <= |value| < B^exp; size == 0
   encodes the value 0 (negative == 0 then).  The number represents a
   ball: the true quantity lies within err ulps of the point value,
   one ulp being B^(exp - size) (for size == 0, B^exp).  Low zero
   limbs are also stripped when err == 0, moving them into exp, so
   exact small integers stay small.

   err is a rigorous bound maintained through every operation, as
   an integer count of ulps of the mantissa's bottom limb.  Products
   and sums compose their bounds in limb arithmetic (128-bit
   intermediate counts at limb anchors, rounded up to whole units);
   divisions and roots go through doubles, rounding up with the fudge
   factor FBALL_EPS.  A bound that reaches B units truncates the
   mantissa by as many limbs as needed, so that err < B and the
   mantissa length tracks the number of accurate limbs; a bound below
   one ulp pads the mantissa with zero limbs down to its own scale.

   Arithmetic functions take a precision n in limbs (the caller
   includes 2-4 guard limbs); mantissas are truncated to about n
   limbs, and less when an operand is accurate to fewer limbs: the
   product of an n/2-accurate and an n/3-accurate ball only ever
   computes ~n/3 limbs. */

typedef struct
{
    nn_ptr d;       /* mantissa limbs */
    slong alloc;    /* allocated limbs */
    slong size;     /* size of mantissa */
    int negative;   /* sign bit */
    slong exp;      /* radix 2^FLINT_BITS exponent */
    ulong err;      /* radius: err ulps of the bottom limb, B^(exp -
                       size) (B^exp for a zero mantissa); 0 = exact.
                       A radius below one ulp is expressed by padding
                       the mantissa with zero limbs down to its scale,
                       one of B ulps or more by truncating the
                       mantissa, so that err < B always holds and at
                       most one limb of noise is carried */
}
fball_struct;

typedef fball_struct fball_t[1];

/* fudge factor: absorbs the rounding of the error bound computations
   that go through doubles (division, roots) */
#define FBALL_EPS (1.0 + 1e-6)

void fball_init(fball_t x);
void fball_clear(fball_t x);
void _fball_grow(fball_t x, slong k);
FLINT_FORCE_INLINE void
fball_fit(fball_t x, slong k)
{
    if (x->alloc < k)
        _fball_grow(x, k);
}

void fball_zero(fball_t x);
void fball_set_ui(fball_t x, ulong c);
void fball_set_si(fball_t x, slong c);
void fball_set(fball_t res, const fball_t x);
void fball_swap(fball_t x, fball_t y);
FLINT_FORCE_INLINE void fball_neg(fball_t x) { if (x->size) x->negative ^= 1; }

int fball_is_zero_exact(const fball_t x);

/* res = a * b, a * c, a + b, a - b at precision n limbs (res may
   alias operands) */
void fball_mul(fball_t res, const fball_t a, const fball_t b, slong n);
void fball_mul_ui(fball_t res, const fball_t a, ulong c, slong n);
/* rr + i ri = (ar + i ai) (br + i bi) to about n limbs in one frame
   (relative to the larger part of the result, both parts kept down to
   the same limb, as the precision of a complex number is that of its
   modulus), by one transform-sharing complex product instead of four;
   outputs may alias inputs */
void fball_mul_complex(fball_t rr, fball_t ri, const fball_t ar, const fball_t ai, const fball_t br, const fball_t bi, slong n);
void fball_div_ui(fball_t res, const fball_t a, ulong c, slong n);
/* res = a - b c to n limbs GIVEN |a - b c| < B^E: only the window of
   the product that reaches the result is computed (a residual of an
   iteration; wrong if the bound fails) */
void fball_submul_bounded(fball_t res, const fball_t a, const fball_t b, const fball_t c, slong E, slong n);
void fball_add(fball_t res, const fball_t a, const fball_t b, slong n);
void fball_sub(fball_t res, const fball_t a, const fball_t b, slong n);
/* res = x + y c resp. x - y c for a word c: in place (res == x) by
   mpn_addmul_1 / mpn_submul_1 in O(|y|) when y c lands inside x's
   window, else y c formed and added */
void fball_addmul_ui(fball_t res, const fball_t x, const fball_t y, ulong c, slong n);
void fball_submul_ui(fball_t res, const fball_t x, const fball_t y, ulong c, slong n);

/* res = a / b resp. 1/sqrt(c) with ~n accurate limbs, by Newton
   iteration (fixed_div_newton / fixed_rsqrt_ui_newton; division by a
   short divisor uses flint_mpn_tdiv_qr instead); b must be nonzero
   with small relative error, 2 <= c < B */
void fball_div(fball_t res, const fball_t a, const fball_t b, slong n);
void fball_rsqrt_ui(fball_t res, ulong c, slong n);
void fball_rsqrt(fball_t res, const fball_t x, slong n);
/* x^(1/k) (k >= 1) and x^(-1/k) (k >= 1), x > 0, k < 2^40 (root.c);
   k = 2 and 3 go to the square and cube roots below, other k to a
   high-order iteration in fball arithmetic whose order r the tuning
   entry point exposes (r = 0: the default) */
void fball_root_ui(fball_t res, const fball_t x, ulong k, slong n);
void fball_rroot_ui(fball_t res, const fball_t x, ulong k, slong n);
void _fball_root_ui_order(fball_t res, const fball_t x, ulong k, slong n, int r, int recip);
void fball_sqrt(fball_t res, const fball_t x, slong n);
void fball_mul_2exp_si(fball_t x, slong e);
/* the arithmetic-geometric mean of x, y >= 0 to about n limbs (agm.c),
   finishing with a series of order m (0 = the default) */
void fball_agm(fball_t res, const fball_t x, const fball_t y, slong n);
void _fball_agm_order(fball_t res, const fball_t x, const fball_t y, slong n, int m);

/* add extra ulps (at the current anchor) to the radius */
void fball_add_error_ulps(fball_t x, double e);
/* add v ulps at the explicit anchor B^anc: use when the intended
   scale is a frame the mantissa may have been normalized away from
   (an exact import strips low zero limbs, so "ulps at the current
   anchor" can silently mean a much coarser frame) */
void fball_add_error(fball_t x, double v, slong anc);
/* enclose also value * t for |t| <= 2^e2 (relative perturbation) */
void fball_add_error_2exp_rel(fball_t x, slong e2);

/* x = (p, len) 2^ebits exactly (p unsigned; len may include zero
   top limbs); the bit part of ebits costs one mpn shift, after
   which all exponent handling is limb-aligned */
void fball_set_mpn_2exp(fball_t x, nn_srcptr p, slong len,
    slong ebits);

/* write a ball known to lie in [0, 1) as a wn-limb fixed-point
   fraction (truncating), returning a rigorous error bound in
   output ulps (2^(-FLINT_BITS wn)); a negative point value within
   the radius is clamped to zero and absorbed into the bound */
double fball_get_fixed(nn_ptr y, slong wn, const fball_t x);

/* write y = floor(x B^n) for a ball known to lie in [0, 1) IF the
   radius determines that floor uniquely, returning 1; returns 0
   otherwise (caller retries at higher precision).  The direct
   verified-floor sibling of fball_get_fixed */
int fball_get_fixed_floor(nn_ptr y, slong n, const fball_t x);

void fball_get_arb(arb_t res, const fball_t x);
void fball_print(const fball_t x);
/* e with |x| < 2^e over the whole ball (value plus radius) */
slong fball_mag_2exp(const fball_t x);
/* e with the relative radius of x (nonzero mantissa) below 2^e */
slong fball_rel_2exp(const fball_t x);
/* x += [-2^e, 2^e] */
void fball_add_error_2exp(fball_t x, slong e);

/* the ball-valued workers of the Newton-Taylor inverses (newton.c):
   -log(x) resp. atan(x) for the fixed-point (x, n); forward and N as
   for the tunable fixed-point entries */
void _fball_neglog_newton(fball_t res, nn_srcptr x, slong n, int forward, slong N);
void _fball_atan_newton(fball_t res, nn_srcptr x, slong n, int forward, slong N);
void _fball_neglog_agm(fball_t res, nn_srcptr x, slong n, slong N);

/* pi to about n limbs (n includes the caller's guard limbs) by
   binary splitting of the Chudnovsky series in fball arithmetic */
void fball_const_pi_chudnovsky(fball_t pi, slong n);

/* atan(p/q) resp. atanh(p/q) (hyperbolic != 0) for 0 <= p < q (p/q at
   most about 0.99) given as mpn integers, to about n limbs (n includes
   the caller's guard limbs), by binary splitting */
void fball_atan_frac_bsplit(fball_t res, nn_srcptr p, slong pn,
    nn_srcptr q, slong qn, int hyperbolic, slong n);

/* Hypergeometric series in the format of y-cruncher's
   SeriesHypergeometric (hypgeom_bsplit.c):

       S = (coefQ + coefP sum_{k>=1} P(k)/Q(k) prod_{j=1}^{k-1} R(j)/Q(j))
           / coefD,

   raised to power = 1 or -1, for integer polynomials P, Q, R
   (coefficient i of degree i; Q(k) != 0 for k >= 1; geometric
   convergence: deg R <= deg Q, |lc R| < |lc Q| if the degrees agree)
   and integers coefP, coefQ, coefD != 0 of any size, given as signed
   mpn integers (-1)^neg (d, n) (n = 0 for zero).  res receives a ball
   for S accurate to about FLINT_BITS (n - 1) bits (n includes the
   caller's guard limbs), by binary splitting in fball arithmetic with a
   rigorous bound for the tail; throws if the tail cannot be bounded.
   The int64 variant takes the coefficients as in y-cruncher's formula
   files. */
typedef struct
{
    nn_srcptr d;
    slong n;
    int neg;
}
fixed_hypgeom_int_struct;

typedef struct
{
    int power;
    fixed_hypgeom_int_struct coefP, coefQ, coefD;
    const fixed_hypgeom_int_struct * P;
    slong Plen;
    const fixed_hypgeom_int_struct * Q;
    slong Qlen;
    const fixed_hypgeom_int_struct * R;
    slong Rlen;
}
fixed_hypgeom_series_struct;

void fball_hypgeom_series(fball_t res, const fixed_hypgeom_series_struct * s,
    slong n);
void fball_hypgeom_series_int64(fball_t res, int power, int64_t coefP,
    int64_t coefQ, int64_t coefD, const int64_t * P, slong Plen,
    const int64_t * Q, slong Qlen, const int64_t * R, slong Rlen, slong n);

/* internal (machin_bsplit.c): res_j = log p_j (the first num primes)
   resp. 2 arg pi_j (the first num nonreal Gaussian primes), j < num,
   to about n limbs; and the exact floors of such values below B into
   the nc-limb entries of a cache (log_primes.c), 0 if undetermined */
void _fixed_log_primes_vec_fball(fball_struct * res, slong num, slong n);
/* internal (machin_bsplit.c): log(u/v), u > v > 0 integers, by Zuniga's
   series [Zun2025] through fball_hypgeom_series */
void _fball_log_ratio_zuniga(fball_t res, const fmpz_t u, const fmpz_t v,
    slong n);
void _fixed_atan_gauss_vec_fball(fball_struct * res, slong num, slong n);
int _fixed_store_floors(nn_ptr e, slong nc, fball_struct * v, slong num);
/* internal thread helpers (parallel.c): run two functions, the second
   on a pool thread if one is free (1 if so), splitting the thread
   budget; run n tasks on up to n threads taking them in order */
int _fixed_parallel_pair(void (* f1)(void *), void * a1,
    void (* f2)(void *), void * a2);
void _fixed_parallel_tasks(void (* f)(slong, void *), void * args, slong n);
/* T = T1 Q2 + P1 T2, Q = Q1 Q2, P = P1 P2 (need_p) in place, T2
   destroyed; par: on two threads */
void _fball_pqt_merge(fball_t P, fball_t Q, fball_t T, fball_t P2,
    fball_t Q2, fball_t T2, int need_p, slong n, int par);
/* the forking rule of the splittings: the halves of a node run on two
   threads (when the budget allows) only if its exact values would stay
   within FIXED_PAR_CAP times the working precision; above, where the
   halves would each hold full-size numbers, they run one after the
   other with the whole budget and the merge runs on two threads */
#define FIXED_PAR_CAP 4.0
void fball_const_log2(fball_t res, slong n);
/* Euler's constant (const_euler.c); set = 0 ... 3 forces log m from
   log 2 / log 2, log 3 / ... / log 2, log 3, log 5, log 7, -1 chooses */
void fball_const_euler(fball_t res, slong n);
/* an fball constant evaluated into an arb ball of prec bits
   (constant_arb.c): the entry point for the arb wrappers */
typedef void (* fixed_constant_func)(fball_t, slong);
void _fixed_constant_arb(arb_t res, fixed_constant_func f, slong prec);

/* the constants of const_e.c, const_log10.c, const_catalan.c,
   const_zeta3.c, const_zeta5.c and const_gamma.c, each computed from
   scratch (no cache) to about FLINT_BITS (n - 1) bits */
void fball_const_e(fball_t res, slong n);
void fball_const_log10(fball_t res, slong n);
void fball_const_catalan(fball_t res, slong n);
void fball_const_zeta3(fball_t res, slong n);
void fball_const_zeta5(fball_t res, slong n);
void fball_const_gamma_1_3(fball_t res, slong n);
void fball_const_gamma_1_4(fball_t res, slong n);
void _fball_const_euler(fball_t res, slong n, int set);

/* verified floor-truncated cached constants: y receives n limbs,
   exactly floor(c B^n) for c = pi/4, log(2) resp. Euler's constant;
   computed limbs are cached per thread and extended on demand (floors nest, so a
   prefix of a longer cached floor is itself the floor) */
void fixed_const_pi_div_4(nn_ptr y, slong n);
void fixed_const_log2(nn_ptr y, slong n);
void fixed_const_euler(nn_ptr y, slong n);
void _fixed_const_pi_div_4_clear(void);
void _fixed_const_log2_clear(void);
void _fixed_const_euler_clear(void);

/* internal: exact-floor entry helpers (tab_exact.c) */
int _fixed_tab_store_floor(nn_ptr e, const arb_t x, slong nc, slong prec);
void _fixed_tab_entry_exact(nn_ptr e, int which, ulong i, slong nc);

/* fast-path entries sit at most a few guard ulps below the truth;
   a guard limb this close to wrapping means the deficit may have
   borrowed into the value limbs, so the entry is recomputed exactly */
#define FIXED_TAB_GUARD_SLACK UWORD(1024)

FLINT_DLL extern const ulong _fixed_exp_logs_static[(FIXED_STATIC_TAB_R + 1) * FIXED_STATIC_TAB_N];
FLINT_DLL extern const ulong _fixed_atans_static[(FIXED_STATIC_TAB_R + 1) * FIXED_STATIC_TAB_N];

/* NOTE: the two inline _tab accessors below dereference the
   thread-local table data and are therefore usable only from
   translation units compiled INTO the library.  Code linking
   against a Windows DLL (tests, tuning and profiling programs)
   must use the exported _fixed_*_entry / _fixed_*_ensure /
   _fixed_*_clear functions instead: thread-local data cannot be
   DLL-exported, so direct references fail to link there. */

FIXED_STATIC_TAB_INLINE nn_srcptr
_fixed_exp_logs_tab(slong nv, slong rc, slong * nc)
{
    if (rc <= FIXED_STATIC_TAB_R && nv + 1 <= FIXED_STATIC_TAB_N)
    {
        *nc = FIXED_STATIC_TAB_N;
        return _fixed_exp_logs_static;
    }
    _fixed_exp_logs_ensure(nv, rc);
    *nc = _fixed_exp_logs_n;
    return _fixed_exp_logs;
}

FIXED_STATIC_TAB_INLINE nn_srcptr
_fixed_atans_tab(slong nv, slong rc, slong * nc)
{
    if (rc <= FIXED_STATIC_TAB_R && nv + 1 <= FIXED_STATIC_TAB_N)
    {
        *nc = FIXED_STATIC_TAB_N;
        return _fixed_atans_static;
    }
    _fixed_atans_ensure(nv, rc);
    *nc = _fixed_atans_n;
    return _fixed_atans;
}

#ifdef __cplusplus
}
#endif

#endif
