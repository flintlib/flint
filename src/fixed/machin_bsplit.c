/*
    Copyright (C) 2013-2014, 2022, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <stdlib.h>
#include "flint.h"
#include "thread_support.h"
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpq.h"
#include "arb.h"
#include "fixed.h"

/* The logarithms of the first num primes and the angles 2 arg(pi_j) of
   the first num nonreal Gaussian primes as fballs, from the
   Machin-type sets of machin_tab.c (the atanh / atan of 1/x_j), one
   followup series each beyond the set, as in arb's
   arb_log_primes_vec_bsplit and arb_atan_gauss_primes_vec_bsplit;
   small precisions read arb's 4608-bit tables of the first 13 values
   directly.

   The atan terms of the Gaussian primes are evaluated by
   fball_atan_frac_bsplit.  The atanh terms of the logarithms are not
   evaluated by their Taylor series but by Zuniga's Ramanujan-type
   series [Zun2025] for log(u/v) = 2 atanh(1/x), u/v = (x+1)/(x-1), in
   the generic fball_hypgeom_series: it gains 6 log2(x) + log2(27/4)
   bits per term against 2 log2(x), while its numbers grow only about
   as much faster, so that the per-term overhead of binary splitting
   (the log2 k bits by which the numbers outgrow the precision gained)
   is about two thirds as large relative to the terms, and much
   smaller for small x (see _zuniga_atanh_inv), from a crossover
   precision growing with x (see _atanh_inv). */

/* ---- Zuniga's series ---- */

/* log(u/v) for integers u > v > 0, by Zuniga's series [Zun2025]

       log(u/v) = sum_{k>=1} rho^k (alpha k + beta) / (gamma k (2k-1))
                  (1)_k (1/2)_k / ((1/6)_k (5/6)_k),

       rho = (u-v)^6 / (108 u^2 v^2 (u+v)^2),
       alpha = -2 (u+v)(u^2 - 14uv + v^2)(u^2 + 4uv + v^2),
       beta = (u+v)^3 (u^2 - 8uv + v^2),  gamma = 2 (u-v)^5,

   (rho = 4 / (27 x^2 (x^2-1)^2) in terms of x), in the format of
   fball_hypgeom_series: with num/den = (u-v)^6 / (6 u^2 v^2 (u+v)^2)
   in lowest terms,

       P = [-(u+v)^2 (u^2 - 8uv + v^2), 2 (u^2 - 14uv + v^2)(u^2 + 4uv + v^2)],
       Q = den [5, -36, 36],  R = num [0, -1, 2],
       coefP / coefD = -(u+v) num / (2 (u-v)^5),  coefQ = 0.

   Used for atanh(1/x) = (1/2) log(u/v), u/v = (x+1)/(x-1) in lowest
   terms (_zuniga_atanh_inv). */
void
_fball_log_ratio_zuniga(fball_t res, const fmpz_t u, const fmpz_t v, slong n)
{
    fmpz_t d, s, t, g, num, den, uu, vv, uv, a, b, zero;
    fmpz P[2], Q[3], R[3];
    fmpq_t c;
    fixed_hypgeom_int_struct ci[11];
    fixed_hypgeom_series_struct ser;
    const fmpz * all[11];
    nn_ptr buf, bp;
    slong i, limbs;

    fmpz_init(d); fmpz_init(s); fmpz_init(t);
    fmpz_init(g); fmpz_init(num); fmpz_init(den); fmpz_init(uu);
    fmpz_init(vv); fmpz_init(uv); fmpz_init(a); fmpz_init(b);
    fmpz_init(zero);
    fmpq_init(c);
    for (i = 0; i < 2; i++)
        fmpz_init(P + i);
    for (i = 0; i < 3; i++)
    {
        fmpz_init(Q + i);
        fmpz_init(R + i);
    }

    fmpz_sub(d, u, v);
    fmpz_add(s, u, v);
    fmpz_mul(uu, u, u);
    fmpz_mul(vv, v, v);
    fmpz_mul(uv, u, v);

    /* num / den */
    fmpz_pow_ui(num, d, 6);
    fmpz_mul(den, uv, s);
    fmpz_mul(den, den, den);
    fmpz_mul_ui(den, den, 6);
    fmpz_gcd(g, num, den);
    fmpz_divexact(num, num, g);
    fmpz_divexact(den, den, g);

    /* P */
    fmpz_add(a, uu, vv);
    fmpz_submul_ui(a, uv, 8);
    fmpz_mul(t, s, s);
    fmpz_mul(P + 0, t, a);
    fmpz_neg(P + 0, P + 0);
    fmpz_add(a, uu, vv);
    fmpz_submul_ui(a, uv, 14);
    fmpz_add(b, uu, vv);
    fmpz_addmul_ui(b, uv, 4);
    fmpz_mul(P + 1, a, b);
    fmpz_mul_2exp(P + 1, P + 1, 1);

    /* Q, R */
    fmpz_mul_ui(Q + 0, den, 5);
    fmpz_mul_si(Q + 1, den, -36);
    fmpz_mul_ui(Q + 2, den, 36);
    fmpz_neg(R + 1, num);
    fmpz_mul_2exp(R + 2, num, 1);

    /* coefP / coefD = -(u+v) num / (2 (u-v)^5) */
    fmpz_mul(t, s, num);
    fmpz_neg(t, t);
    fmpz_pow_ui(g, d, 5);
    fmpz_mul_2exp(g, g, 1);
    fmpq_set_fmpz_frac(c, t, g);

    /* as signed mpn integers */
    all[0] = P; all[1] = P + 1;
    all[2] = Q; all[3] = Q + 1; all[4] = Q + 2;
    all[5] = R; all[6] = R + 1; all[7] = R + 2;
    all[8] = fmpq_numref(c); all[9] = zero; all[10] = fmpq_denref(c);
    limbs = 0;
    for (i = 0; i < 11; i++)
        limbs += FLINT_MAX(fmpz_size(all[i]), 1);
    bp = buf = flint_malloc(limbs * sizeof(ulong));
    for (i = 0; i < 11; i++)
    {
        slong m = fmpz_size(all[i]);
        fmpz_abs(t, all[i]);
        if (m > 0)
            fmpz_get_ui_array(bp, m, t);
        ci[i].d = bp;
        ci[i].n = m;
        ci[i].neg = (fmpz_sgn(all[i]) < 0);
        bp += FLINT_MAX(m, 1);
    }

    ser.power = 1;
    ser.P = ci;
    ser.Plen = 2;
    ser.Q = ci + 2;
    ser.Qlen = 3;
    ser.R = ci + 5;
    ser.Rlen = 3;
    ser.coefP = ci[8];
    ser.coefQ = ci[9];
    ser.coefD = ci[10];

    fball_hypgeom_series(res, &ser, n);

    flint_free(buf);
    fmpz_clear(d); fmpz_clear(s); fmpz_clear(t);
    fmpz_clear(g); fmpz_clear(num); fmpz_clear(den); fmpz_clear(uu);
    fmpz_clear(vv); fmpz_clear(uv); fmpz_clear(a); fmpz_clear(b);
    fmpz_clear(zero);
    fmpq_clear(c);
    for (i = 0; i < 2; i++)
        fmpz_clear(P + i);
    for (i = 0; i < 3; i++)
    {
        fmpz_clear(Q + i);
        fmpz_clear(R + i);
    }
}

/* atanh(1/x) = (1/2) log(u/v), u/v = (x+1)/(x-1) in lowest terms */
static void
_zuniga_atanh_inv(fball_t res, nn_srcptr x, slong xn, slong n)
{
    fmpz_t u, v;

    fmpz_init(u);
    fmpz_init(v);
    fmpz_set_ui_array(u, x, xn);
    fmpz_sub_ui(v, u, 1);
    fmpz_add_ui(u, u, 1);
    if (fmpz_is_even(u))
    {
        fmpz_fdiv_q_2exp(u, u, 1);
        fmpz_fdiv_q_2exp(v, v, 1);
    }
    _fball_log_ratio_zuniga(res, u, v, n);
    fball_mul_2exp_si(res, -1);
    fmpz_clear(u);
    fmpz_clear(v);
}

/* atanh(1/x): Zuniga's series from about n >= bits(x)^2 / 2 limbs, the
   Taylor series below (measured crossovers: x ~ 2^12 always, 2^30 at
   about 400 limbs, 2^64 at 2000, 2^77 at 3000; the larger x, the
   longer the terms of Zuniga's series and the more precision it takes
   for the saving on the per-term overhead to pay for them) */
static void
_atanh_inv(fball_t res, nn_srcptr x, slong xn, slong n)
{
    ulong one = 1;
    slong b = (xn - 1) * FLINT_BITS + FLINT_BIT_COUNT(x[xn - 1]);

    if (2 * n >= b * b)
        _zuniga_atanh_inv(res, x, xn, n);
    else
        fball_atan_frac_bsplit(res, &one, 1, x, xn, 1, n);
}

/* ---- the static tables ---- */

#define STATIC_NUM 13
#define STATIC_LIMBS (4608 / FLINT_BITS)

/* x = table entry (a fraction of STATIC_LIMBS limbs) times 2^e, within
   one ulp of the table */
static void
_set_static(fball_t x, nn_srcptr tab, slong e)
{
    slong t = e - FLINT_BITS * STATIC_LIMBS;
    slong q = t >> (FLINT_BITS == 64 ? 6 : 5);

    fball_set_mpn_2exp(x, tab, STATIC_LIMBS, t);
    fball_add_error(x, ldexp(1.0, (int) (t - q * FLINT_BITS)), q);
}

/* ---- the combination with the Machin coefficients ---- */

/* acc (an limbs, room for one more) += (x, xn) (c, cn) */
static slong
_acc_addmul(nn_ptr acc, slong an, nn_srcptr x, slong xn, nn_srcptr c,
    slong cn, nn_ptr tmp)
{
    ulong cy;
    slong tn;

    if (cn == 1)
    {
        cy = mpn_addmul_1(acc, x, xn, c[0]);
        cy = mpn_add_1(acc + xn, acc + xn, an - xn, cy);
    }
    else
    {
        flint_mpn_mul(tmp, x, xn, c, cn);
        tn = xn + cn;
        cy = mpn_add(acc, acc, an, tmp, tn);
    }
    FLINT_ASSERT(cy == 0);
    (void) cy;
    return an;
}

/* res_i = (1/den) sum_j c_ij y_j for i < min(num, tab->num), the y_j
   in [0, 1).  The y_j are read once as wn-limb fractions Y_j with
   error bounds e_j (fball_get_fixed); the dot products are exact mpn
   sums (nonnegative and negative coefficients accumulated apart), with
   radius sum |c_ij| e_j; one division by den follows.  The y_j are
   freed (left zero). */
static void
_machin_combine(fball_struct * res, slong num, fball_struct * y,
    const fixed_machin_struct * tab, slong wn)
{
    slong i, j, ln = tab->num, climbs, an;
    nn_ptr Y, pos, neg, tmp;
    double * ey;
    fball_t t, d;

    fball_init(t);
    fball_init(d);
    fball_set_ui(d, tab->den);

    climbs = (tab->cbits + FLINT_BITS) / FLINT_BITS + 1;
    an = wn + climbs + FLINT_BIT_COUNT(ln) / FLINT_BITS + 2;

    Y = flint_malloc((ln * wn + 3 * an) * sizeof(ulong));
    pos = Y + ln * wn;
    neg = pos + an;
    tmp = neg + an;
    ey = flint_malloc(ln * sizeof(double));

    /* the y_j are released as they are read (only the fixed-point
       copies are needed from here) */
    for (j = 0; j < ln; j++)
    {
        ey[j] = fball_get_fixed(Y + j * wn, wn, y + j);
        fball_clear(y + j);
        fball_init(y + j);
    }

    for (i = 0; i < FLINT_MIN(num, ln); i++)
    {
        const unsigned char * csign;
        nn_srcptr crow = fixed_machin_c_row_raw(tab, i, &climbs, &csign);
        double err = 0.0;

        flint_mpn_zero(pos, an);
        flint_mpn_zero(neg, an);

        for (j = 0; j < ln; j++)
        {
            nn_srcptr c = crow + j * climbs;
            slong cn = climbs;

            while (cn > 0 && c[cn - 1] == 0)
                cn--;
            if (cn == 0)
                continue;
            _acc_addmul(csign[j] ? neg : pos, an, Y + j * wn, wn, c, cn,
                tmp);
            /* |c| < (top + 1) B^(cn - 1) */
            err += ldexp((double) c[cn - 1] + 1.0, FLINT_BITS * (cn - 1))
                * (1.0 + 0x1p-50) * ey[j];
        }

        /* res = (pos - neg) B^-wn +/- err ulps */
        if (mpn_cmp(pos, neg, an) >= 0)
        {
            mpn_sub_n(pos, pos, neg, an);
            fball_set_mpn_2exp(res + i, pos, an, -FLINT_BITS * wn);
        }
        else
        {
            mpn_sub_n(pos, neg, pos, an);
            fball_set_mpn_2exp(res + i, pos, an, -FLINT_BITS * wn);
            fball_neg(res + i);
        }
        fball_add_error(res + i, err * (1.0 + 1e-10) + 1.0, -wn);

        if (tab->den != 1)
        {
            fball_div(t, res + i, d, wn);
            fball_swap(res + i, t);
        }
    }

    flint_free(Y);
    flint_free(ey);
    fball_clear(t);
    fball_clear(d);
}

/* one series of a set: atanh(1/x) (kind 0), atan(1/x) (kind 1) or
   atan(p/x) (kind 2) */
typedef struct
{
    fball_struct * res;
    nn_srcptr x;
    slong xn;
    ulong p;
    int kind;
    double cost;
}
machin_job_t;

typedef struct
{
    machin_job_t * jobs;
    slong n;
}
machin_work_t;

static void
_machin_worker(slong i, void * arg)
{
    machin_work_t * w = (machin_work_t *) arg;
    const machin_job_t * J = w->jobs + i;
    ulong one = 1;

    if (J->kind == 0)
        _atanh_inv(J->res, J->x, J->xn, w->n);
    else
        fball_atan_frac_bsplit(J->res, J->kind == 2 ? &J->p : &one, 1,
            J->x, J->xn, 0, w->n);
}

static void
_machin_job(machin_job_t * J, fball_struct * res, nn_srcptr x, slong xn,
    ulong p, int kind)
{
    double lx;

    while (xn > 1 && x[xn - 1] == 0)
        xn--;
    J->res = res;
    J->x = x;
    J->xn = xn;
    J->p = p;
    J->kind = kind;
    /* terms, about n / log2(x/p) */
    lx = (xn - 1) * (double) FLINT_BITS + log2((double) x[xn - 1] + 1.0);
    if (kind == 2)
        lx -= log2((double) p);
    J->cost = 1.0 / FLINT_MAX(lx, 0.01);
}

static int
_machin_job_cmp(const void * a, const void * b)
{
    double ca = ((const machin_job_t *) a)->cost;
    double cb = ((const machin_job_t *) b)->cost;
    return (ca < cb) - (ca > cb);
}

/* all the series of a set at n limbs, on the available threads, taken
   costliest first: at most one series per thread is in flight, and a
   thread left over when there are fewer series goes to the splitting of
   one of them */
static void
_machin_run(machin_job_t * jobs, slong num, slong n)
{
    machin_work_t w;

    qsort(jobs, num, sizeof(machin_job_t), _machin_job_cmp);
    w.jobs = jobs;
    w.n = n;
    _fixed_parallel_tasks(_machin_worker, &w, num);
}

static fball_struct *
_machin_y_init(slong num)
{
    slong j;
    fball_struct * y = flint_malloc(num * sizeof(fball_struct));
    for (j = 0; j < num; j++)
        fball_init(y + j);
    return y;
}

static void
_machin_series_clear(fball_struct * y, slong num)
{
    slong j;
    for (j = 0; j < num; j++)
        fball_clear(y + j);
    flint_free(y);
}

/* guard limbs for the combination: |c_ij| < 2^(cbits - 1) amplifies
   the absolute errors of the y_j (and the terms cancel) */
static slong
_machin_guard(const fixed_machin_struct * tab)
{
    return (tab->cbits + FLINT_BIT_COUNT(tab->num) + FLINT_BITS - 1)
        / FLINT_BITS + 1;
}

void
_fixed_log_primes_vec_fball(fball_struct * res, slong num, slong n)
{
    const fixed_machin_struct * tab;
    fball_struct * y;
    slong i, j, k, ln, wp;
    ulong * primes, * qv;
    n_primes_t iter;
    fball_t t;
    machin_job_t * jobs;

    FLINT_ASSERT(num >= 1 && num <= FIXED_LOG_PRIMES_MAX);

    if (num <= STATIC_NUM && n + 1 < STATIC_LIMBS)
    {
        for (i = 0; i < num; i++)
            _set_static(res + i, arb_log_p_tab[i],
                (i >= 1) + (i >= 4) + (i >= 16));
        return;
    }

    tab = fixed_machin_table(0, num);
    ln = tab->num;
    wp = n + _machin_guard(tab);

    /* the series: atanh(1/x_j) of the table and, for the followups
       log p = atanh(1/(2p^2 - 1)) + (1/2) (log((p-1)/2) + log((p+1)/2))
       + log 2 (the factors of (p -/+ 1)/2 being smaller primes), the
       atanh(1/(2p^2 - 1)) directly into res */
    primes = flint_malloc(2 * num * sizeof(ulong));
    qv = primes + num;
    n_primes_init(iter);
    for (i = 0; i < num; i++)
        primes[i] = n_primes_next(iter);
    n_primes_clear(iter);

    y = _machin_y_init(ln);
    jobs = flint_malloc(FLINT_MAX(num, ln) * sizeof(machin_job_t));
    for (j = 0; j < ln; j++)
        _machin_job(jobs + j, y + j, tab->x + j * FIXED_MACHIN_X_LIMBS,
            FIXED_MACHIN_X_LIMBS, 1, 0);
    for (i = ln; i < num; i++)
    {
        qv[i] = 2 * primes[i] * primes[i] - 1;
        _machin_job(jobs + ln + (i - ln), res + i, qv + i, 1, 1, 0);
    }
    _machin_run(jobs, ln + FLINT_MAX(num - ln, 0), wp);
    flint_free(jobs);

    _machin_combine(res, num, y, tab, wp);
    _machin_series_clear(y, ln);

    fball_init(t);
    for (i = ln; i < num; i++)
    {
        ulong p = primes[i], h;
        int side;

        fball_mul_2exp_si(res + i, 1);

        for (side = 0; side < 2; side++)
        {
            n_factor_t fac;

            h = side ? (p + 1) / 2 : (p - 1) / 2;
            n_factor_init(&fac);
            n_factor(&fac, h, 1);
            for (j = 0; j < fac.num; j++)
                for (k = 0; k < i; k++)
                    if (fac.p[j] == primes[k])
                    {
                        fball_addmul_ui(res + i, res + i, res + k, fac.exp[j], wp);
                    }
        }

        fball_mul_2exp_si(res + i, -1);
        fball_add(res + i, res + i, res + 0, wp);
    }
    fball_clear(t);
    flint_free(primes);
}

void
_fixed_atan_gauss_vec_fball(fball_struct * res, slong num, slong n)
{
    const fixed_machin_struct * tab;
    fball_struct * y;
    slong i, j, ln, wp;
    slong * fol;
    machin_job_t * jobs;

    FLINT_ASSERT(num >= 1 && num <= FIXED_ATAN_GAUSS_MAX);

    if (num <= STATIC_NUM && n + 1 < STATIC_LIMBS)
    {
        /* 2 arg(pi_i) */
        static const char exponents[STATIC_NUM] =
            { 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1 };
        for (i = 0; i < num; i++)
            _set_static(res + i, arb_atan_gauss_tab[i], exponents[i] + 1);
        return;
    }

    tab = fixed_machin_table(1, num);
    ln = tab->num;
    wp = n + _machin_guard(tab);

    /* the series: atan(1/x_j) of the table and, for the followups
       arg pi_i = atan(t) + arg pi_j with the neighbour j (among the
       previous ones) minimizing |t|, t = (xb ya - xa yb) / (xa ya +
       xb yb), atan(|t|) directly into res */
    y = _machin_y_init(ln);
    jobs = flint_malloc(FLINT_MAX(num, ln) * sizeof(machin_job_t));
    fol = flint_malloc(3 * FLINT_MAX(num - ln, 1) * sizeof(slong));
    for (j = 0; j < ln; j++)
        _machin_job(jobs + j, y + j, tab->x + j * FIXED_MACHIN_X_LIMBS,
            FIXED_MACHIN_X_LIMBS, 1, 1);
    for (i = ln; i < num; i++)
    {
        double best = 100, tt;
        slong xa, xb, ya, yb, best_j = 0, p, q;
        slong * f = fol + 3 * (i - ln);

        xa = _fixed_gaussian_primes[2 * i];
        xb = _fixed_gaussian_primes[2 * i + 1];

        for (j = 0; j < FLINT_MIN(i, 100); j++)
        {
            ya = _fixed_gaussian_primes[2 * j];
            yb = _fixed_gaussian_primes[2 * j + 1];
            tt = (xb * ya - xa * yb) / (double) (xa * ya + xb * yb);
            if (fabs(tt) < best)
            {
                best = fabs(tt);
                best_j = j;
            }
        }

        ya = _fixed_gaussian_primes[2 * best_j];
        yb = _fixed_gaussian_primes[2 * best_j + 1];
        p = xb * ya - xa * yb;
        q = xa * ya + xb * yb;
        FLINT_ASSERT(q > 0 && p != 0);
        f[0] = best_j;
        f[1] = p;
        f[2] = q;
        _machin_job(jobs + ln + (i - ln), res + i, (nn_srcptr) (f + 2), 1,
            FLINT_UABS(p), 2);
    }
    _machin_run(jobs, ln + FLINT_MAX(num - ln, 0), wp);
    flint_free(jobs);

    _machin_combine(res, num, y, tab, wp);
    _machin_series_clear(y, ln);

    for (i = ln; i < num; i++)
    {
        slong * f = fol + 3 * (i - ln);
        if (f[1] < 0)
            fball_neg(res + i);
        fball_add(res + i, res + i, res + f[0], wp);
    }
    flint_free(fol);

    /* theta = 2 arg */
    for (i = 0; i < num; i++)
        fball_mul_2exp_si(res + i, 1);
}
