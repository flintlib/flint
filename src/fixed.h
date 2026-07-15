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
   4 = full bit-burst.  Error at most FIXED_EXP_REDUCED_MAX_ERR
   ulps. */
#define FIXED_EXP_REDUCED_MAX_ERR 96
void fixed_exp_reduced(nn_ptr y, nn_srcptr t, slong wn,
    flint_bitcnt_t r, int alg);

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
   values by mpn binary splitting (at most the true value, short by a
   couple of ulps); i >= 1. */
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
   |error| <= 4 B^-n / sqrt(a). */
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
void _fixed_exp_sum_bs_powtab(nn_ptr T, slong * tn, nn_ptr Q,
    slong * qn, flint_bitcnt_t * Qexp, nn_srcptr xp, slong xn,
    flint_bitcnt_t r, slong N);

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
