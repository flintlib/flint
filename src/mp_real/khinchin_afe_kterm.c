/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "flint.h"
#include "longlong.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "mag.h"
#include "mp_real.h"
#include "impl.h"

/* The k-terms of the theta splitting of Khinchin's constant
   (const_khinchin_afe.c): with K = k^2, x = K u0, u0 = 4^-a,

       I_k = int_(u0)^inf e^(-K u) E(4u) / u du
           = e^(-x) sum_(n>=1) (4/K)^n e_n(x) / n^2,

   e_n(x) = sum_(i<n) x^i / i!.  With P_n = (4/K)^n e_n(x),
   Q_n = (4 u0)^n / n! and U_n = P_n / n,

       U_(n+1) = (4/K) (n U_n + Q_n) / (n + 1),
       Q_(n+1) = 4 u0 Q_n / (n + 1),     I_k = e^(-x) sum_n U_n / n,

   a triangular 2 x 2 recurrence with one sum over the common
   denominator q_n = K (n + 1) n 4^a, summed by the matrix binary
   splitting of _mp_real_linrec (linrec.c): a node over [a, b) is the
   transition [A B; 0 D] with the sum row (Sp, Sq) over Q, and two nodes
   combine as one augmented product

       [A B; 0 D; Sp Sq] = [A2 B2 0; 0 D2 0; Sp2 Sq2 Q2] [A1 B1; 0 D1; Sp1 Sq1],

   which shares the transforms of the operands among the products
   (measured faster than the separate products with an mpz content
   removal of the weights 1/n, a form of the content removal of Cheng,
   Hanrot, Thome, Zima and Zimmermann, which saved 5% at most).  The
   exact products of chunks of about half the working precision are
   applied to the state in ball arithmetic.  Subtrees of at most KM_LEAVES
   leaves are done in fixed-size limb arrays on the stack (all entries
   are nonnegative), without fmpz or mpz overhead (a fifth of the time
   in fmpz).

   TRUNCATION.  Term n of e^-x sum P_n / n^2 is
   t_n = (4/K)^n Pr[Poisson(x) <= n - 1] / n^2, bounded by
   b_n = (4/K)^n min(1, Pi_n) / n^2 with Pi_n = n x^(n-1) e^-x / (n-1)!
   for n - 1 < x (the Poisson weights increase up to x) and by
   (4/K)^n / n^2 for n - 1 >= x.  In the first regime the ratio of
   consecutive bounds is at most 4 u0 / (n + 1) <= 1/2 (u0 <= 1/4), in
   the second at most 4/K <= 4/9, so with
   n0 = ceil(x) + 1 the tail from N is at most
   2 b_N + (9/5) (4/K)^n0 / n0^2 (N < n0) resp. (9/5) b_N (N >= n0). */

/* limb-array bottom ********************************************************/

/* A leaf entry is a product of at most four words below 2^40 (K),
   2^26 (n, n + 1) and 2^16 (4 u0d), at most 108 bits, and a node over l
   leaves at most l times the bits of a leaf plus a carry per level. */
#define KM_LEAVES 32
#define KM_LIMBS ((128 * KM_LEAVES) / FLINT_BITS + 4)
#define KM_N_MAX (WORD(1) << 26)
#if FLINT_BITS == 64
#define KM_K_MAX (UWORD(1) << 40)
#else
#define KM_K_MAX UWORD_MAX
#endif
#define KM_U0D_MAX (UWORD(1) << 14)

typedef struct
{
    ulong e[6][KM_LIMBS];
    slong n[6];
}
kmnode_struct;

enum { KM_A = 0, KM_B, KM_D, KM_SP, KM_SQ, KM_Q };

/* r = a b c d, normalized (a zero factor gives 0) */
static void
_km_set4(ulong * r, slong * rn, ulong a, ulong b, ulong c, ulong d)
{
    slong n = 1;

    if (a == 0 || b == 0 || c == 0 || d == 0)
    {
        *rn = 0;
        return;
    }
    r[0] = a;
    r[n] = mpn_mul_1(r, r, n, b); n += (r[n] != 0);
    r[n] = mpn_mul_1(r, r, n, c); n += (r[n] != 0);
    r[n] = mpn_mul_1(r, r, n, d); n += (r[n] != 0);
    *rn = n;
}

/* r = a b (normalized; r must not alias a or b) */
static slong
_km_mul(nn_ptr r, nn_srcptr a, slong an, nn_srcptr b, slong bn)
{
    slong rn;

    if (an == 0 || bn == 0)
        return 0;
    if (an >= bn)
        flint_mpn_mul(r, a, an, b, bn);
    else
        flint_mpn_mul(r, b, bn, a, an);
    rn = an + bn;
    while (rn > 0 && r[rn - 1] == 0)
        rn--;
    return rn;
}

/* r += t (r has room for max(rn, tn) + 1 limbs) */
static slong
_km_add(nn_ptr r, slong rn, nn_srcptr t, slong tn)
{
    ulong cy;

    if (tn == 0)
        return rn;
    if (rn == 0)
    {
        flint_mpn_copyi(r, t, tn);
        return tn;
    }
    if (rn >= tn)
    {
        cy = mpn_add(r, r, rn, t, tn);
        r[rn] = cy;
        return rn + (cy != 0);
    }
    cy = mpn_add(r, t, tn, r, rn);
    r[tn] = cy;
    return tn + (cy != 0);
}

static void
_km_get_fmpz(fmpz_t x, nn_srcptr e, slong n)
{
    if (n == 0)
        fmpz_zero(x);
    else
        fmpz_set_ui_array(x, e, n);
}

static void
km_leaf(kmnode_struct * N, slong n, ulong K, ulong u0d)
{
    _km_set4(N->e[KM_Q], &N->n[KM_Q], K, n + 1, n, u0d);
    _km_set4(N->e[KM_B], &N->n[KM_B], 4, n, u0d, 1);
    _km_set4(N->e[KM_A], &N->n[KM_A], 4 * u0d, n, n, 1);
    _km_set4(N->e[KM_D], &N->n[KM_D], 4, K, n, 1);
    _km_set4(N->e[KM_SP], &N->n[KM_SP], K, n + 1, u0d, 1);
    N->n[KM_SQ] = 0;
}

static void
km_combine(kmnode_struct * R, const kmnode_struct * N1, const kmnode_struct * N2)
{
    ulong t[KM_LIMBS];
    slong tn;
#define E1(i) N1->e[i], N1->n[i]
#define E2(i) N2->e[i], N2->n[i]
    /* Sq = Q2 Sq1 + Sp2 B1 + Sq2 D1, Sp = Q2 Sp1 + Sp2 A1 */
    R->n[KM_SQ] = _km_mul(R->e[KM_SQ], E2(KM_Q), E1(KM_SQ));
    tn = _km_mul(t, E2(KM_SP), E1(KM_B));
    R->n[KM_SQ] = _km_add(R->e[KM_SQ], R->n[KM_SQ], t, tn);
    tn = _km_mul(t, E2(KM_SQ), E1(KM_D));
    R->n[KM_SQ] = _km_add(R->e[KM_SQ], R->n[KM_SQ], t, tn);
    R->n[KM_SP] = _km_mul(R->e[KM_SP], E2(KM_Q), E1(KM_SP));
    tn = _km_mul(t, E2(KM_SP), E1(KM_A));
    R->n[KM_SP] = _km_add(R->e[KM_SP], R->n[KM_SP], t, tn);
    /* B = A2 B1 + B2 D1, A = A1 A2, D = D1 D2, Q = Q1 Q2 */
    R->n[KM_B] = _km_mul(R->e[KM_B], E2(KM_A), E1(KM_B));
    tn = _km_mul(t, E2(KM_B), E1(KM_D));
    R->n[KM_B] = _km_add(R->e[KM_B], R->n[KM_B], t, tn);
    R->n[KM_A] = _km_mul(R->e[KM_A], E1(KM_A), E2(KM_A));
    R->n[KM_D] = _km_mul(R->e[KM_D], E1(KM_D), E2(KM_D));
    R->n[KM_Q] = _km_mul(R->e[KM_Q], E1(KM_Q), E2(KM_Q));
#undef E1
#undef E2
}

static void
km_bsplit(kmnode_struct * R, slong a, slong b, ulong K, ulong u0d)
{
    if (b - a == 1)
    {
        km_leaf(R, a, K, u0d);
    }
    else
    {
        kmnode_struct N1, N2;
        slong m = a + (b - a) / 2;
        km_bsplit(&N1, a, m, K, u0d);
        km_bsplit(&N2, m, b, K, u0d);
        km_combine(R, &N1, &N2);
    }
}

/* the recurrence for _mp_real_linrec ****************************************/

typedef struct
{
    ulong K, u0d;
}
_kterm_ctx;

/* state (U, Q), sum row (1/n, 0): the matrix [A B; 0 D; Sp 0] over q */
static void
_kterm_leaf(fmpz_mat_t T, fmpz_t phi, fmpz_t q, slong n, const void * vctx)
{
    const _kterm_ctx * c = (const _kterm_ctx *) vctx;
    fmpz * Sp = fmpz_mat_entry(T, 2, 0);

    fmpz_set_ui(Sp, c->K);
    fmpz_mul_ui(Sp, Sp, n + 1);
    fmpz_mul_ui(Sp, Sp, c->u0d);
    fmpz_mul_ui(q, Sp, n);
    fmpz_set_ui(fmpz_mat_entry(T, 0, 1), 4);
    fmpz_mul_ui(fmpz_mat_entry(T, 0, 1), fmpz_mat_entry(T, 0, 1), n);
    fmpz_mul_ui(fmpz_mat_entry(T, 0, 1), fmpz_mat_entry(T, 0, 1), c->u0d);
    fmpz_mul_ui(fmpz_mat_entry(T, 0, 0), fmpz_mat_entry(T, 0, 1), n);
    fmpz_set_ui(fmpz_mat_entry(T, 1, 1), c->K);
    fmpz_mul_ui(fmpz_mat_entry(T, 1, 1), fmpz_mat_entry(T, 1, 1), 4);
    fmpz_mul_ui(fmpz_mat_entry(T, 1, 1), fmpz_mat_entry(T, 1, 1), n);
}

static void
_kterm_bottom(fmpz_mat_t T, fmpz_t phi, fmpz_t q, slong a, slong b, const void * vctx)
{
    const _kterm_ctx * c = (const _kterm_ctx *) vctx;
    kmnode_struct M;

    km_bsplit(&M, a, b, c->K, c->u0d);
    _km_get_fmpz(fmpz_mat_entry(T, 0, 0), M.e[KM_A], M.n[KM_A]);
    _km_get_fmpz(fmpz_mat_entry(T, 0, 1), M.e[KM_B], M.n[KM_B]);
    _km_get_fmpz(fmpz_mat_entry(T, 1, 1), M.e[KM_D], M.n[KM_D]);
    _km_get_fmpz(fmpz_mat_entry(T, 2, 0), M.e[KM_SP], M.n[KM_SP]);
    _km_get_fmpz(fmpz_mat_entry(T, 2, 1), M.e[KM_SQ], M.n[KM_SQ]);
    _km_get_fmpz(q, M.e[KM_Q], M.n[KM_Q]);
}

/* the number of terms for a tail below 2^-(prec + 16) (doubles: an
   estimate only; the tail actually added is bounded rigorously) */
static slong
_kterm_terms(ulong K, double u0, slong prec)
{
    double x = (double) K * u0, n0 = ceil(x) + 1, second, tb, lq = log2(4.0 / (double) K);
    slong N;

    second = log2(1.8) + n0 * lq - 2 * log2(n0);
    for (N = FLINT_MAX(2, (slong) (8 * u0) + 2); ; N++)
    {
        tb = N * lq - 2 * log2((double) N);
        if (N - 1 < x)
        {
            tb += FLINT_MIN(0.0, (-x + (N - 1) * log(x) - lgamma((double) N)
                + log((double) N)) * 1.4426950408889634);
            tb = FLINT_MAX(1 + tb, second) + 1;
        }
        else
            tb += log2(1.8);
        if (tb < -(double) prec - 16)
            return N;
        /* while the second part is too large, no N < n0 works */
        if (N - 1 < x && second >= -(double) prec - 16 && N < (slong) n0 - 1)
            N = (slong) n0 - 1;
    }
}

/* the tail bound (see the top of the file) */
static void
_kterm_tail(mag_t res, ulong K, slong N, slong a)
{
    mag_t q, b, s, x, t;
    ulong n0;

    mag_init(q); mag_init(b); mag_init(s); mag_init(x); mag_init(t);

    /* x = K u0 (rounded down), n0 = ceil(x) + 1 */
    mag_set_ui_lower(x, K);
    mag_mul_2exp_si(x, x, -2 * a);
    n0 = (K >> (2 * a)) + ((K & ((UWORD(1) << (2 * a)) - 1)) != 0) + 1;

    mag_set_ui(q, 4);
    mag_div_ui(q, q, K);
    mag_pow_ui(b, q, N);
    mag_div_ui(b, b, N);
    mag_div_ui(b, b, N);

    if ((ulong) N < n0)     /* N - 1 < x */
    {
        /* Pi_N = N x^(N-1) e^-x / (N-1)! */
        mag_t pi;
        mag_init(pi);
        mag_set_ui(t, K);
        mag_mul_2exp_si(t, t, -2 * a);
        mag_pow_ui(pi, t, N - 1);
        mag_expinv(t, x);
        mag_mul(pi, pi, t);
        mag_mul_ui(pi, pi, N);
        mag_rfac_ui(t, N - 1);
        mag_mul(pi, pi, t);
        if (mag_cmp_2exp_si(pi, 0) < 0)
            mag_mul(b, b, pi);
        mag_mul_2exp_si(b, b, 1);
        /* + (9/5) (4/K)^n0 / n0^2 */
        mag_pow_ui(s, q, n0);
        mag_div_ui(s, s, n0);
        mag_div_ui(s, s, n0);
        mag_mul_ui(s, s, 9);
        mag_div_ui(s, s, 5);
        mag_add(b, b, s);
        mag_clear(pi);
    }
    else
    {
        mag_mul_ui(b, b, 9);
        mag_div_ui(b, b, 5);
    }
    mag_swap(res, b);

    mag_clear(q); mag_clear(b); mag_clear(s); mag_clear(x); mag_clear(t);
}

/* the exact levels up to this many times the working precision (measured
   at 3 10^5 digits: 1, 2, 4, 8 within 10%, 4 best) */

void
_mp_real_khinchin_afe_kterm(mp_real_t res, ulong k,
    const _mp_real_khinchin_afe_struct * ctx, const mp_real_t ek)
{
    ulong K = k * k;
    double u0 = ldexp(1.0, -2 * ctx->a);
    slong prec = ctx->prec, wp, N, n;
    _kterm_ctx kc;
    _mp_real_linrec_struct rec;
    mp_real_struct x0[2];
    mp_real_t t;
    mag_t tail;

    if (k < 3 || k >= (UWORD(1) << (FLINT_BITS / 2)))
        flint_throw(FLINT_ERROR, "_mp_real_khinchin_afe_kterm: k out of range\n");

    mag_init(tail);

    /* I_k is about e^-((K - 4) u0): the relative precision needed */
    wp = prec - (slong) ((K - 4.0) * u0 * 1.4426950408889634) + 32;
    if (wp < 32)
    {
        /* 0 < I_k <= e^(-(K-4) u0) / ((K-4) u0), as E(4u)/u <= e^(4u)/u0
           for u >= u0 */
        mag_t x;
        mag_init(x);
        mag_set_ui_lower(x, K - 4);
        mag_mul_2exp_si(x, x, -2 * ctx->a);
        mag_expinv(tail, x);
        mag_div_ui(tail, tail, K - 4);
        mag_mul_2exp_si(tail, tail, 2 * ctx->a);
        mp_real_zero(res);
        _mp_real_add_error_mag(res, tail);
        mag_clear(x);
        mag_clear(tail);
        return;
    }

    N = _kterm_terms(K, u0, prec);
    n = mp_real_prec_bits(wp + 32);

    kc.K = K;
    kc.u0d = ctx->u0d;
    rec.d = 2;
    rec.s = 1;
    rec.g = 0;
    rec.leaf = _kterm_leaf;
    rec.ctx = &kc;
    if (N < KM_N_MAX && K < KM_K_MAX && ctx->u0d <= KM_U0D_MAX)
    {
        rec.bottom = _kterm_bottom;
        rec.bottom_len = KM_LEAVES;
    }
    else
    {
        rec.bottom = NULL;
        rec.bottom_len = 0;
    }

    mp_real_init(x0 + 0);
    mp_real_init(x0 + 1);
    mp_real_init(t);
    mp_real_set_ui(x0 + 0, 4);
    mp_real_div_ui(x0 + 0, x0 + 0, K, n);
    mp_real_set_ui(x0 + 1, 1);
    mp_real_mul_2exp_si(x0 + 1, x0 + 1, 2 - 2 * ctx->a);
    /* sum_(n >= 1) U_n / n from (U_1, Q_1) = (4/K, 4 u0) */
    mp_real_zero(res);
    _mp_real_linrec_eval(x0, res, NULL, &rec, 1, N, NULL, 1, mp_real_prec_bits(wp + 64));

    if (ek != NULL)
        mp_real_mul(res, res, ek, mp_real_prec_bits(wp));
    else
    {
        mp_real_set_ui(t, K);
        mp_real_mul_2exp_si(t, t, -2 * ctx->a);
        mp_real_neg(t, t);
        mp_real_exp_bits(t, t, wp + 32);
        mp_real_mul(res, res, t, mp_real_prec_bits(wp));
    }

    _kterm_tail(tail, K, N, ctx->a);
    _mp_real_add_error_mag(res, tail);

    mp_real_clear(x0 + 0);
    mp_real_clear(x0 + 1);
    mp_real_clear(t);
    mag_clear(tail);
}
