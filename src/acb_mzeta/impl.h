/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Internal declarations of the acb_mzeta module: the module's translation
   units include this as "impl.h", the tests and profiling programs as
   "acb_mzeta/impl.h". */

#ifndef ACB_MZETA_IMPL_H
#define ACB_MZETA_IMPL_H

#include <math.h>
#include <stdint.h>
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_mat.h"
#include "arb.h"
#include "acb.h"
#include "arb_mat.h"
#include "acb_mzeta.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Cache ********************************************************************/

/* entry types */
#define _ACB_MZETA_CACHE_LATTICE 1
#define _ACB_MZETA_CACHE_FTAB 2

/* the most recent entries of each type are kept */
#define _ACB_MZETA_CACHE_MAX_PER_TYPE 8

typedef int (* _acb_mzeta_cache_match_func)(const void * obj, const void * key);
typedef void (* _acb_mzeta_cache_clear_func)(void * obj);

void * _acb_mzeta_cache_find(acb_mzeta_cache_t cache, int type,
    _acb_mzeta_cache_match_func match, const void * key);
void _acb_mzeta_cache_insert(acb_mzeta_cache_t cache, int type, void * obj,
    _acb_mzeta_cache_clear_func clear);

/* Reduced lattice ***********************************************************/

/* The lattice A Z^d after LLL reduction of the Gram matrix, with the data
   of the splitting at lam = |det A|^(1/d) shared by the Crandall-type
   algorithms: the quadratic forms of both sums as Cholesky factors,
   t = |Cd (n - u)|^2 = (pi/lam^2) |A (n - u)|^2 and
   t = |Cr (m + w)|^2 = pi lam^2 |A^-T (m + w)|^2, and the constants of
   the theta tail bound (acb_theta_sum_radius) for both. */
typedef struct
{
    slong d;
    slong prec;
    int ok;
    arb_mat_t A0;           /* the input basis (cache key) */
    arb_mat_t A, Ainv;      /* reduced basis */
    arb_mat_t G, Ginv;      /* Gram matrix A^T A and its inverse */
    arb_t vol, lam, lam2;
    arb_mat_t Cd, Cr;       /* upper triangular */
    arb_t Kd, Kr;           /* (1 + sqrt(8/pi)) prod_j (1 + sqrt(2 pi) / C_jj) */
}
_acb_mzeta_lattice_struct;

/* upper triangular C with C^T C = M; the constant
   (1 + sqrt(8/pi)) prod_j (1 + sqrt(2 pi)/C_jj) of the theta tail bound */
int _acb_mzeta_upper_cho(arb_mat_t C, const arb_mat_t M, slong prec);
void _acb_mzeta_theta_tail_const(arb_t K, const arb_mat_t C, slong prec);

const _acb_mzeta_lattice_struct * _acb_mzeta_lattice_get(acb_mzeta_cache_t cache,
    const arb_mat_t A, slong prec);

/* Epstein zeta: the reduced problem *****************************************/

/* x = A (u0 + m), A^T y = w = w0 + l with m, l integer and u0, w0 near
   zero; Z(x; y) = e(-w.m) Z(A u0; A^-T w0), e(t) = exp(2 pi i t).
   The direct sum runs over n with t = |Cd (n - u0)|^2 and phase e(-w0.n),
   omitting n = 0 if x is in the lattice; the reciprocal sum over m with
   t = |Cr (m + w0)|^2 and phase e(-u0.(m + w0)), omitting m = 0 if y is in
   the dual lattice and, for the regularised function, m = l (the term of
   y itself, replaced by -g_{d-nu}). */
typedef struct
{
    slong d;
    arb_ptr u0, w0, w;
    fmpz * m, * l;
    int xs, ys;             /* x in L, y in L* */
    int regflag;            /* the regularised function was requested */
    int reg;                /* ... and the term m = l is replaced (y != 0) */
    slong * lpt;            /* l as slongs (reg only) */
    acb_t nu, a, b;         /* nu, nu/2, (d - nu)/2 */
    int real;               /* nu is real */
    double v0, v1;          /* max(0, Re a - 1), max(0, Re b - 1) (+ slack) */
}
_acb_mzeta_epstein_problem_struct;

int _acb_mzeta_epstein_problem_init(_acb_mzeta_epstein_problem_struct * P,
    const _acb_mzeta_lattice_struct * L, const acb_t nu, arb_srcptr x, arb_srcptr y,
    int flags, slong prec);
void _acb_mzeta_epstein_problem_clear(_acb_mzeta_epstein_problem_struct * P);

/* whether the point n of the direct (side 0) or reciprocal (side 1) sum is
   omitted */
int _acb_mzeta_epstein_skip(const _acb_mzeta_epstein_problem_struct * P,
    int side, const slong * n);

/* truncation: the smallest T >= max(4, v + 1) with
   K T^((d-1)/2) e^-T / (T - v) <= 2^-bits, and a rigorous bound for it
   (the theta tail bound times |G| <= e^-t / (t - v)) */
double _acb_mzeta_epstein_choose_T(const arb_t K, slong d, double v, slong bits);
void _acb_mzeta_epstein_tail_bound(mag_t res, const arb_t K, slong d, double T, double v);

/* res = Z from the direct and reciprocal sums (balls including the
   truncation errors): the regularisation term, the singular terms, the
   prefactor (pi/lam^2)^(nu/2)/Gamma(nu/2) and the translation phase */
void _acb_mzeta_epstein_finish(acb_t res, const acb_t Sd, const acb_t Sr,
    const _acb_mzeta_lattice_struct * L, const _acb_mzeta_epstein_problem_struct * P,
    slong prec, slong wp);

/* Incomplete gamma ***********************************************************/

/* F(t) = e^t t^-a Gamma(a, t) by the asymptotic series, if it reaches the
   relative tolerance 2^-prec; returns 0 otherwise (res untouched).  t exact,
   t > 0. */
int _acb_mzeta_gamma_upper_scaled_asymp(acb_t res, const acb_t a, const arb_t t, slong prec);

/* Tables of F(t) = e^t t^-a Gamma(a, t) for the double-precision kernels:
   Taylor expansions of degree K on intervals |t - t0| <= rho with centers
   decreasing geometrically from tmax to tmin, coefficients as n-term
   double expansions (n = 1, 2), with bounds E (value) and D (derivative):
   for |t - t0| + r <= rho, |F(t') - p(t - t0)| <= E + D r for |t' - t| <= r,
   p the Horner value computed by _acb_mzeta_ftab_eval. */
typedef struct
{
    double a;               /* real index */
    int n;                  /* double expansion length */
    slong K;                /* degree */
    slong K1;               /* n = 2: degrees >= K1 in plain double */
    slong len;              /* number of intervals */
    double * t0, * rho, * E, * D;
    double * c;             /* len (K + 1) n coefficients */
    double tlo, thi;        /* the covered range */
    slong * lut;            /* first interval per bucket of t */
    slong kmin, kmax;
    double tmin_req, tmax_req;
}
_acb_mzeta_ftab_struct;

void _acb_mzeta_ftab_init(_acb_mzeta_ftab_struct * T, double a, int n, double tmin, double tmax);
void _acb_mzeta_ftab_clear(_acb_mzeta_ftab_struct * T);
const _acb_mzeta_ftab_struct * _acb_mzeta_ftab_get(acb_mzeta_cache_t cache, double a, int n,
    double tmin, double tmax);

/* F at len points: t, res are arrays of dNb balls (n doubles and a
   radius); entries whose t is outside the table get radius +inf.
   The Horner evaluation runs in plain dN arithmetic with an a priori
   error bound. */
void _acb_mzeta_ftab_eval(double * res, const double * t, slong len,
    const _acb_mzeta_ftab_struct * T);

/* Shells of a shifted integral lattice ***************************************/

/* For an integral positive definite F (d x d, d <= 8), sigma, ell in Q^d
   (denominators <= 64): V(n) = F(S n - s), s = S sigma, takes the values
   g M + r0 (M >= 0); the counts are C_r(M) = #{n in Z^d : V(n) = g M + r0,
   N ell.n = r (mod N)} for M < len.

   They are produced by handlers (shells_*.c), each of which recognises
   the cases it can count and estimates its cost; compute takes the
   cheapest applicable handler (falling back to the next one if it fails)
   and the consumers only see the counts block by block.  Costs are
   estimates of the work in nanoseconds, independent of the number of
   threads, so that the choices (and hence the results) are too. */

#define ACB_MZETA_SHELLS_MAX_DIM 8
#define ACB_MZETA_SHELLS_NUM_HANDLERS 4

typedef struct _acb_mzeta_shells_struct _acb_mzeta_shells_struct;

typedef struct
{
    const char * name;
    /* handler data if the handler applies to T (after the common part of
       T is set), else NULL */
    void * (*init)(const _acb_mzeta_shells_struct * T);
    void (*clear)(void * data);
    /* estimated work for len shells, INFINITY if not possible */
    double (*cost)(const void * data, const _acb_mzeta_shells_struct * T, ulong len);
    /* prepares the counts of the shells M < len; 0 on failure */
    int (*compute)(void * data, const _acb_mzeta_shells_struct * T, ulong len);
    /* out[i N + r] = C_r(M0 + i), 0 <= i < n; thread-safe */
    void (*block)(ulong * out, const void * data, const _acb_mzeta_shells_struct * T,
        ulong M0, ulong n);
    /* an upper bound for the counts */
    ulong (*max_count)(const void * data, const _acb_mzeta_shells_struct * T);
    /* mask[r] = 0 if C_r(M) = 0 for all M < len (optional: NULL = all
       residues may occur) */
    void (*residues)(char * mask, const void * data, const _acb_mzeta_shells_struct * T);
}
_acb_mzeta_shells_handler_struct;

struct _acb_mzeta_shells_struct
{
    slong d;
    fmpz_mat_t F;
    fmpz * s, * ell;        /* s = S sigma; ell = N ell mod N */
    ulong S, N, g, r0;
    ulong len;
    int method;             /* 0 automatic, else 1 + the forced handler */
    int variant;            /* handler options for testing (0 = default) */
    slong active;           /* the handler of the last compute, or -1 */
    void * data[ACB_MZETA_SHELLS_NUM_HANDLERS];
};

extern const _acb_mzeta_shells_handler_struct _acb_mzeta_shells_r3;
extern const _acb_mzeta_shells_handler_struct _acb_mzeta_shells_r3_odd;
extern const _acb_mzeta_shells_handler_struct _acb_mzeta_shells_products;
extern const _acb_mzeta_shells_handler_struct _acb_mzeta_shells_enum;

/* in this order (the method numbers) */
#define ACB_MZETA_SHELLS_R3 1
#define ACB_MZETA_SHELLS_R3_ODD 2
#define ACB_MZETA_SHELLS_PRODUCTS 3
#define ACB_MZETA_SHELLS_ENUM 4

/* init returns 0 (and leaves nothing to clear) when the denominators or
   the values are too large */
int _acb_mzeta_shells_init(_acb_mzeta_shells_struct * T, const fmpz_mat_t F,
    const fmpq * sigma, const fmpq * ell);
void _acb_mzeta_shells_clear(_acb_mzeta_shells_struct * T);
double _acb_mzeta_shells_cost(const _acb_mzeta_shells_struct * T, ulong len);
int _acb_mzeta_shells_compute(_acb_mzeta_shells_struct * T, ulong len);
void _acb_mzeta_shells_block(ulong * out, const _acb_mzeta_shells_struct * T, ulong M0, ulong n);
ulong _acb_mzeta_shells_max_count(const _acb_mzeta_shells_struct * T);
void _acb_mzeta_shells_residues(char * mask, const _acb_mzeta_shells_struct * T);
const char * _acb_mzeta_shells_name(const _acb_mzeta_shells_struct * T);

/* F(x) = sum_i q_ii (x_i + sum_{j>i} q_ij x_j)^2 in double (for bounds
   and estimates only); returns 0 if not positive */
int _acb_mzeta_ldl_d(double * q, const fmpz_mat_t F);

/* helpers for handlers that store all counts C[M N + r], M < len */
void _acb_mzeta_shells_table_block(ulong * out, const ulong * C, ulong N, ulong M0, ulong n);
ulong _acb_mzeta_shells_table_max(const ulong * C, ulong len, ulong N);
void _acb_mzeta_shells_table_residues(char * mask, const ulong * C, ulong len, ulong N);

/* Multiple zeta values *******************************************************/

/* the binary word x0^(s_1 - 1) x1 ... x0^(s_r - 1) x1 (0 = x0, 1 = x1) of a
   composition; returns its length (the weight) */
slong _acb_mzeta_mzv_word(unsigned char * w, const ulong * s, slong r);

/* the iterated integral of an admissible word (w[0] = 0, w[k-1] = 1) by
   Akhilesh's double tails in fixed point */
void _acb_mzeta_mzv_akhilesh(arb_t res, const unsigned char * w, slong k, slong prec);

/* the same by Hoelder convolution and binary splitting (quasilinear) */
void _acb_mzeta_mzv_bsplit(arb_t res, const unsigned char * w, slong k, slong prec);

/* Goncharov's G(a_l, ..., a_{L-1}; z) for l = 0, ..., L (res[L] = 1),
   a_{L-1} != 0, |z| < 0.95 min |a_l| over the nonzero letters. Uses the
   Taylor series at zero; returns 0 if not applicable. */
int _acb_mzeta_G_suffixes(acb_ptr res, acb_srcptr a, slong L, const acb_t z,
    slong prec);

/* G(a; 1) directly or by Hoelder convolution at 1/2; 0 if not applicable */
int _acb_mzeta_G_one(acb_t res, acb_srcptr a, slong L, slong prec);

/* the same for nw admissible words at once, sharing the subword states */
void _acb_mzeta_mzv_akhilesh_vec(arb_ptr res, const unsigned char * const * words,
    const slong * lens, slong nw, slong prec);

/* Madelung: shared building blocks ******************************************/

/* Deterministic shared accumulator: every contribution is rounded to the
   grid 2^-E on its own and added as an exact integer, with the radius in
   whole ulps, so the sum does not depend on the order of the additions
   (the number of threads) and only one full-precision number is stored. */
typedef struct
{
    void * lock;
    fmpz_t mid, rad;
    slong E;
}
_acb_mzeta_fixacc_struct;

void _acb_mzeta_fixacc_init(_acb_mzeta_fixacc_struct * A, slong E);
void _acb_mzeta_fixacc_clear(_acb_mzeta_fixacc_struct * A);
void _acb_mzeta_fixacc_add(_acb_mzeta_fixacc_struct * A, const arb_t x, int negate);
void _acb_mzeta_fixacc_get(arb_t res, const _acb_mzeta_fixacc_struct * A);

/* r3(N) for N < X; r3(3 + 8k) for k < K (one FFT each); t2(j) for j < K
   (segmented divisor sieve, in parallel); r3(3 + 8k), k0 <= k < k0 + len,
   from t2 by sparse convolution */
uint32_t * _acb_mzeta_r3_all_table(ulong X);   /* r3(N), N < X */
uint32_t * _acb_mzeta_r3_3mod8_table(ulong K);
uint16_t * _acb_mzeta_t2_table(ulong K, slong nthreads);
void _acb_mzeta_r3_3mod8_block(uint32_t * out, const uint16_t * t2, ulong k0, ulong len);

/* Fixed-point sums with small rational coefficients: on (fn + 3)-limb
   fixed-point numbers (scaled by B^fn), bs = sum_i q_i x_i +
   sum_chunks floor(sum_{i in chunk} r_i (D/den_i) x_i / D), where
   num_i = q_i den_i + r_i, x_i is the (fn + 1)-limb entry at
   xtab + i stride, and the chunks are formed so that r_i D/den_i fits in
   limbs (1 or 2) limbs; the divisions of several chunks are interleaved.
   Requires den_i >= 1.  Given an error of at most xerr ulps (B^-fn) in
   each x_i, returns a bound for the error of bs in ulps:
   sum_i (q_i + 1) xerr + the number of chunks (integer arithmetic). */
ulong _acb_mzeta_fixed_ratsum(ulong * bs, const ulong * xtab, slong stride, slong fn,
    const ulong * num, const ulong * den, ulong len, ulong xerr, int limbs,
    void * chunks, ulong * work);
/* the scratch sizes for _acb_mzeta_fixed_ratsum with len terms */
slong _acb_mzeta_fixed_ratsum_chunks_size(ulong len);
slong _acb_mzeta_fixed_ratsum_work_size(slong fn);

/* P such sums over the same x_i at once (sum p: num + p len, den + p len,
   limbs[p], result at bs + p (fn + 3), error bound err[p]), going through
   the terms once so that each x_i is read from memory once; chunks has P
   times the chunks size */
void _acb_mzeta_fixed_ratsum_multi(ulong * bs, ulong * err, const ulong * xtab, slong stride,
    slong fn, const ulong * num, const ulong * den, ulong len, slong P, ulong xerr,
    const int * limbs, void * chunks, ulong * work);
slong _acb_mzeta_fixed_ratsum_multi_work_size(slong fn, slong P);
/* 1 or 2: the cheaper divisor size for these terms */
int _acb_mzeta_fixed_ratsum_choose_limbs(const ulong * num, const ulong * den, ulong len);

#ifdef __cplusplus
}
#endif

#endif
