/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "ulong_extras.h"
#include "fmpz_vec.h"
#include "fmpq_vec.h"
#include "fmpz_mat.h"
#include "fmpq_mat.h"
#include "arb_hypgeom.h"
#include "thread_support.h"
#include "impl.h"

/* the precision below which the irrational side uses the series rather
   than bit-burst steps (measured) */
#define SERIES_MAX_TARGET 16384

/* The Epstein zeta function of an integral lattice at rational points by
   shells (the generalization of the Madelung algorithm).

   For an integral positive definite Gram matrix F and u, v in Q^d,

     Z = sum'_n e(-v.n) F(n - u)^(-s),  s = nu/2,

   the Crandall splitting with a free parameter kappa gives

     Z = (1/Gamma(s)) [kappa^s S_D
            + pi^(d/2) det(F)^(-1/2) kappa^(s - d/2) e(-u.v) S_R]  + singular terms,
     S_D = sum'_n e(-v.n) G_s(kappa F(n - u)),
     S_R = sum'_m e(-u.m) G_b(pi^2 F^-1(m + v) / kappa),  b = d/2 - s,

   G_a(t) = t^-a Gamma(a, t); singular terms -kappa^s / Gamma(s + 1) if
   u is in Z^d and pi^(d/2) det(F)^(-1/2) kappa^(s - d/2) / ((s - d/2) Gamma(s))
   if v is in Z^d.  Both sums are sums over shells (_acb_mzeta_shells with
   the forms F and F_R = delta F^-1, the shifts u and -v, the characters v
   and u), with exponents t = gamma V, V = g M + r0 the shell values.

   kappa is chosen so that one side (B) has t = lambda V with a dyadic
   lambda: its incomplete gamma values come from binary splitting over
   exact numbers.  The other side (O) has t = pi^2 V / (lambda Pi),
   Pi = S_D^2 S_R^2 delta.  When its index is a positive integer, it is a
   short sum of power series in q = e^(-gamma g) with rational
   coefficients C / V^j (as in the Madelung code: a shared fixed-point
   table of powers of q, denominators cleared by _acb_mzeta_fixed_ratsum);
   otherwise its incomplete gamma values come from the series of
   _acb_mzeta_gamma_upper_scaled or, at high precision, a bit-burst step
   each (gamma_bitburst), and lambda puts few shells there.
   lambda balances the shell counts: B costs ~T / (lambda g_B) incomplete
   gamma evaluations, O ~T lambda Pi / (pi^2 g_O) series terms.

   Truncation: each side is cut at a shell value whose tail is bounded by
   the theta tail bound (_acb_mzeta_epstein_tail_bound) for the Cholesky
   factor of gamma S^2 F; |G_a(t)| <= e^-t / (t - max(0, a - 1)). */

#define SIDE_QSERIES 0
#define SIDE_BSPLIT 1
#define SIDE_GENERIC 2

#define MIN_PREC 30
#define TERM_GUARD 8

typedef struct
{
    int kind;
    _acb_mzeta_shells_struct sh;
    fmpq_t a;               /* index */
    slong sint;             /* QSERIES: a as an integer */
    ulong len;              /* shells M < len */
    arb_t gam;              /* t = gam V (lambda for BSPLIT) */
    double gamd;
    int skip0;              /* omit V = 0 (r0 = 0, M = 0) */
    slong cbits;            /* bits of the largest count */
    acb_ptr chi;            /* e(-r/N) */
    acb_t sum;
    mag_t tail;
    slong target;           /* absolute accuracy 2^-target */
}
side_struct;

static slong
taper(slong prec, double log2_mag)
{
    slong p = prec + TERM_GUARD + (slong) ceil(log2_mag);
    return FLINT_MAX(p, MIN_PREC);
}

/* weight of a shell with the counts C[r]: sum_r C_r chi_r */
static void
shell_weight(acb_t w, const side_struct * Sd, const ulong * C, slong prec)
{
    ulong N = Sd->sh.N, r;
    acb_t t;

    if (N == 1)
    {
        acb_set_ui(w, C[0]);
        return;
    }
    if (N == 2)
    {
        acb_set_si(w, (slong) C[0] - (slong) C[1]);
        return;
    }
    acb_init(t);
    acb_zero(w);
    for (r = 0; r < N; r++)
        if (C[r] != 0)
        {
            acb_mul_ui(t, Sd->chi + r, C[r], prec);
            acb_add(w, w, t, prec);
        }
    acb_clear(t);
}

static ulong
shell_total(const side_struct * Sd, const ulong * C)
{
    ulong N = Sd->sh.N, r, s = 0;
    for (r = 0; r < N; r++)
        s += C[r];
    return s;
}

/* Incomplete gamma side (BSPLIT and GENERIC) *********************************/

typedef struct
{
    side_struct * Sd;
    arb_srcptr eml;         /* e^-lambda (BSPLIT) */
    slong nblocks;
    _acb_mzeta_fixacc_struct * acc_re, * acc_im;
    int series;             /* GENERIC: the convergent/asymptotic series */
    acb_srcptr aa, ga;      /* the index and Gamma of it (series) */
}
gamma_ctx;

/* e^-z to relative precision p: eml^V, the given value e0, or exp */
static void
expmz_get(arb_t e, const arb_t z, const arb_t eml, ulong V, const arb_t e0, slong p)
{
    if (e0 != NULL)
        arb_set_round(e, e0, p);
    else if (eml != NULL)
    {
        arb_set_round(e, eml, p + 2 * FLINT_BIT_COUNT(V) + 8);
        arb_pow_ui(e, e, V, p + 2 * FLINT_BIT_COUNT(V) + 8);
    }
    else
    {
        arb_neg(e, z);
        arb_exp(e, e, p + 8);
    }
}

/* G_a(z) for exact z > 0, with e^-z = eml^V (eml != NULL), e0 (given to
   sufficient precision, see gamma_bitburst) or computed; |error| <=
   2^-tolbits */
static void
gamma_exact(arb_t G, const fmpq_t a, const arb_t z, const arb_t eml, ulong V,
    const arb_t e0, slong tolbits)
{
    mag_t tol, err, za;
    arb_t t, u, expmz;
    slong N, p;
    double zd = arf_get_d(arb_midref(z), ARF_RND_DOWN), ad = fmpq_get_d(a);

    mag_init(tol);
    mag_init(err);
    mag_init(za);
    arb_init(t);
    arb_init(u);
    arb_init(expmz);

    /* tolerance for Gamma(a, z): 2^-tolbits z^a (z^a bounded below) */
    {
        fmpq_t ma;
        fmpq_init(ma);
        arb_pow_fmpq(t, z, a, 30);
        arb_get_mag_lower(tol, t);
        mag_mul_2exp_si(tol, tol, -tolbits);
        fmpq_neg(ma, a);
        arb_pow_fmpq(t, z, ma, 30);
        arb_get_mag(za, t);    /* z^-a */
        fmpq_clear(ma);
    }

    N = _arb_hypgeom_gamma_upper_fmpq_inf_choose_N(err, a, z, tol);
    if (N != -1)
    {
        /* Gamma(a, z) = z^(a-1) e^-z S, so G = e^-z S / z */
        p = taper(tolbits, -zd * 1.4426950408889634 - log2(zd));
        expmz_get(expmz, z, eml, V, e0, p);
        _arb_hypgeom_gamma_upper_fmpq_inf_bsplit_sum(G, a, z, N, p);
        arb_mul(G, G, expmz, p);
        arb_div(G, G, z, p);
        mag_mul(err, err, za);
        arb_add_error_mag(G, err);
    }
    else if (!(fmpz_is_one(fmpq_denref(a)) && fmpz_sgn(fmpq_numref(a)) <= 0))
    {
        /* gamma(a, z) = z^a e^-z S / a: G = z^-a Gamma(a) - e^-z S / a,
           both terms of size z^-a Gamma(a) */
        double lg = (-ad * log(zd) + lgamma(fabs(ad) + (ad <= 0 ? 1.0 : 0.0))) * 1.4426950408889634;
        p = tolbits + 16 + (slong) FLINT_MAX(0.0, lg + 4 + fabs(ad) * 2);
        N = _arb_hypgeom_gamma_lower_fmpq_0_choose_N(err, a, z, tol);
        expmz_get(expmz, z, eml, V, e0, p);
        _arb_hypgeom_gamma_lower_fmpq_0_bsplit_sum(G, a, z, N, p);
        arb_mul(G, G, expmz, p);
        arb_div_fmpz(G, G, fmpq_numref(a), p);
        arb_mul_fmpz(G, G, fmpq_denref(a), p);
        arb_gamma_fmpq(t, a, p);
        {
            fmpq_t ma;
            fmpq_init(ma);
            fmpq_neg(ma, a);
            arb_pow_fmpq(u, z, ma, p);
            fmpq_clear(ma);
        }
        arb_mul(t, t, u, p);
        arb_sub(G, t, G, p);
        mag_mul(err, err, za);
        arb_add_error_mag(G, err);
    }
    else
    {
        /* a = -n: Gamma(-n, z) by the convergent series (cancellation
           ~ e^z), G = z^n Gamma(-n, z) */
        slong n = -fmpz_get_si(fmpq_numref(a));
        mag_t tol2;
        mag_init(tol2);
        arb_pow_ui(t, z, n, 30);
        arb_get_mag_lower(tol2, t);
        mag_one(err);
        mag_div(tol2, err, tol2);
        mag_mul_2exp_si(tol2, tol2, -tolbits);
        p = tolbits + 16 + (slong) (zd * 2 * 1.4426950408889634 + n * log2(zd + 1.0));
        N = _arb_hypgeom_gamma_upper_singular_si_choose_N(err, n, z, tol2);
        arb_set_round(u, z, p);
        _arb_hypgeom_gamma_upper_singular_si_bsplit(G, n, u, N, p);
        arb_add_error_mag(G, err);
        arb_pow_ui(t, z, n, p);
        arb_mul(G, G, t, p);
        mag_clear(tol2);
    }

    mag_clear(tol);
    mag_clear(err);
    mag_clear(za);
    arb_clear(t);
    arb_clear(u);
    arb_clear(expmz);
}

/* G_a(z) for an inexact z > 0 (bit-burst): G at a dyadic z0 of b bits
   near z by gamma_exact, then the Taylor series in eps = z - z0, with the
   coefficients g_k = G^(k)(z0)/k! from the differential equation
   t G'' + (1 + a + t) G' + a G = 0, g_1 = -(e^-z0 + a G(z0)) / z0:

     g_{k+2} = -((k+1)(k+1+a+z0) g_{k+1} + (k+a) g_k) / (z0 (k+1)(k+2)).

   G(t) = int_1^oo v^(a-1) e^(-t v) dv gives |g_k| <= Gamma(k+a) / (k! z0^(k+a))
   for k + a > 0, so the tail beyond K is bounded by its first term over
   1 - q, q = max(1, (K+a)/(K+1)) |eps| / z0. */
static void
gamma_bitburst(arb_t G, const fmpq_t a, const arb_t z, slong tolbits)
{
    slong b, K, k, p;
    double zd = arf_get_d(arb_midref(z), ARF_RND_DOWN), ad = fmpq_get_d(a), lg;
    arb_t z0, eps, e0, g0, g1, g2, t, ak;
    mag_t epsm, z0l, term, q, tail;

    b = FLINT_MAX(64, (slong) (2.0 * sqrt((double) tolbits)));

    arb_init(z0); arb_init(eps); arb_init(e0); arb_init(g0); arb_init(g1); arb_init(g2);
    arb_init(t); arb_init(ak);
    mag_init(epsm); mag_init(z0l); mag_init(term); mag_init(q); mag_init(tail);

    arf_set_round(arb_midref(z0), arb_midref(z), b, ARF_RND_DOWN);
    mag_zero(arb_radref(z0));
    arb_sub(eps, z, z0, tolbits + 64);
    arb_get_mag(epsm, eps);
    arb_get_mag_lower(z0l, z0);

    /* e^-z0 to the precision gamma_exact needs in either branch */
    lg = (-ad * log(zd) + lgamma(fabs(ad) + 1.0)) * 1.4426950408889634;
    p = tolbits + 64 + (slong) FLINT_MAX(0.0, lg + 4 + fabs(ad) * 2);
    arb_neg(e0, z0);
    arb_exp(e0, e0, p);
    gamma_exact(g0, a, z0, NULL, 0, e0, tolbits + 8);

    /* K and the tail bound: term_k = Gamma(k+a)/k! |eps|^k / z0^(k+a) */
    for (K = 1; ; K++)
    {
        if (K + ad > 0.5)
        {
            fmpq_t ka;
            fmpq_init(ka);
            fmpq_add_si(ka, a, K);
            arb_gamma_fmpq(t, ka, 30);
            arb_get_mag(term, t);
            mag_rfac_ui(q, K);
            mag_mul(term, term, q);
            mag_pow_ui(q, epsm, K);
            mag_mul(term, term, q);
            arb_pow_fmpq(t, z0, ka, 30);
            arb_get_mag_lower(q, t);
            mag_div(term, term, q);
            /* q = max(1, (K+a)/(K+1)) |eps| / z0 */
            mag_div(q, epsm, z0l);
            if (ad > 1)
            {
                mag_mul_ui(q, q, K + (slong) ceil(ad));
                mag_div_ui(q, q, K + 1);
            }
            fmpq_clear(ka);
            if (mag_cmp_2exp_si(q, -1) < 0)
            {
                mag_geom_series(tail, q, 0);
                mag_mul(tail, tail, term);
                if (mag_cmp_2exp_si(tail, -tolbits - 4) < 0)
                    break;
            }
        }
        if (K > 100000)
        {
            mag_inf(tail);
            break;
        }
    }

    /* coefficients and Horner (relative precision for the size of G) */
    p = taper(tolbits, -zd * 1.4426950408889634 - log2(zd + 1.0)) + FLINT_BIT_COUNT(K) + 8;
    {
        arb_ptr g = _arb_vec_init(K);
        arb_set_round(g + 0, g0, p);
        if (K > 1)
        {
            /* g_1 = -(e^-z0 + a g_0) / z0 */
            arb_set_fmpq(ak, a, p);
            arb_mul(g1, ak, g + 0, p);
            arb_add(g1, g1, e0, p);
            arb_div(g1, g1, z0, p);
            arb_neg(g + 1, g1);
        }
        for (k = 0; k + 2 < K; k++)
        {
            /* (k+1)(k+1+a+z0) g_{k+1} + (k+a) g_k */
            arb_add_ui(t, ak, k + 1, p);
            arb_add(t, t, z0, p);
            arb_mul_ui(t, t, k + 1, p);
            arb_mul(g2, t, g + k + 1, p);
            arb_add_ui(t, ak, k, p);
            arb_addmul(g2, t, g + k, p);
            arb_div(g2, g2, z0, p);
            arb_div_ui(g2, g2, (k + 1) * (k + 2), p);
            arb_neg(g + k + 2, g2);
        }
        arb_zero(G);
        for (k = K - 1; k >= 0; k--)
        {
            arb_mul(G, G, eps, p);
            arb_add(G, G, g + k, p);
        }
        _arb_vec_clear(g, K);
    }
    arb_add_error_mag(G, tail);

    arb_clear(z0); arb_clear(eps); arb_clear(e0); arb_clear(g0); arb_clear(g1); arb_clear(g2);
    arb_clear(t); arb_clear(ak);
    mag_clear(epsm); mag_clear(z0l); mag_clear(term); mag_clear(q); mag_clear(tail);
}

static void
gamma_worker(slong b, void * vctx)
{
    gamma_ctx * C = vctx;
    side_struct * Sd = C->Sd;
    arb_t z, e, G;
    acb_t w, s;
    ulong M;

    arb_init(z);
    arb_init(e);
    arb_init(G);
    acb_init(w);
    acb_init(s);

    for (M = b; M < Sd->len; M += C->nblocks)
    {
        ulong cnt[64];
        ulong V = Sd->sh.g * M + Sd->sh.r0, tot;
        slong p, tb;
        double t;

        _acb_mzeta_shells_block(cnt, &Sd->sh, M, 1);
        tot = shell_total(Sd, cnt);

        if (tot == 0 || (V == 0 && Sd->skip0))
            continue;

        t = Sd->gamd * V;
        /* each of the len terms to 2^-(target + 4) / len */
        tb = Sd->target + FLINT_BIT_COUNT(tot) + FLINT_BIT_COUNT(Sd->len) + 4;
        p = taper(tb, -t * 1.4426950408889634 - log2(t));

        if (Sd->kind == SIDE_BSPLIT)
        {
            arb_set_ui(z, V);
            arb_mul(z, z, Sd->gam, ARF_PREC_EXACT);
            gamma_exact(G, Sd->a, z, C->eml, V, NULL, tb);
            acb_set_arb(s, G);
        }
        else
        {
            arb_mul_ui(z, Sd->gam, V, tb + 64);
            if (C->series)
                _acb_mzeta_gamma_upper_scaled(s, C->aa, z, C->ga, FLINT_MAX(p, 32) + 8);
            else
            {
                gamma_bitburst(G, Sd->a, z, tb);
                acb_set_arb(s, G);
            }
        }

        shell_weight(w, Sd, cnt, p);
        acb_mul(s, s, w, p);
        _acb_mzeta_fixacc_add(C->acc_re, acb_realref(s), 0);
        _acb_mzeta_fixacc_add(C->acc_im, acb_imagref(s), 0);
    }

    arb_clear(z);
    arb_clear(e);
    arb_clear(G);
    acb_clear(w);
    acb_clear(s);
}

static void
gamma_side(side_struct * Sd, const arb_t lam, slong nthreads)
{
    gamma_ctx C;
    _acb_mzeta_fixacc_struct re, im;
    arb_t eml;
    acb_t aa, ga;

    arb_init(eml);
    acb_init(aa);
    acb_init(ga);

    /* irrational arguments: below SERIES_MAX_TARGET bits the series of
       _acb_mzeta_gamma_upper_scaled (rectangular splitting) is faster than
       a bit-burst step (measured 1.3-2x for the whole sum); above, the
       bit-burst step, whose cost grows more slowly with the precision */
    C.series = (Sd->kind != SIDE_BSPLIT && Sd->target < SERIES_MAX_TARGET);
    if (C.series)
    {
        /* Gamma(a) to about the precision the series needs (it is
           recomputed there otherwise) */
        arb_set_fmpq(acb_realref(aa), Sd->a, Sd->target + 256);
        acb_gamma(ga, aa, Sd->target + 3 * FLINT_BITS + 64);
    }
    C.aa = aa;
    C.ga = ga;
    if (Sd->kind == SIDE_BSPLIT)
    {
        arb_neg(eml, lam);
        arb_exp(eml, eml, Sd->target + 128);
    }
    _acb_mzeta_fixacc_init(&re, Sd->target + 32);
    _acb_mzeta_fixacc_init(&im, Sd->target + 32);

    C.Sd = Sd;
    C.eml = eml;
    C.nblocks = FLINT_MAX(1, FLINT_MIN(1024, (slong) Sd->len));
    C.acc_re = &re;
    C.acc_im = &im;
    flint_parallel_do(gamma_worker, &C, C.nblocks, nthreads, FLINT_PARALLEL_STRIDED);

    _acb_mzeta_fixacc_get(acb_realref(Sd->sum), &re);
    _acb_mzeta_fixacc_get(acb_imagref(Sd->sum), &im);
    _acb_mzeta_fixacc_clear(&re);
    _acb_mzeta_fixacc_clear(&im);
    arb_clear(eml);
    acb_clear(aa);
    acb_clear(ga);
}

/* Power series side ***********************************************************/

typedef struct
{
    side_struct * Sd;
    slong nres;
    ulong * res;            /* residues with nonzero counts */
    const ulong * ftab;     /* q^i, i <= m, entry i at ftab + i (FN + 1) */
    slong FN;
    ulong ftab_err;
    arb_srcptr qm;
    ulong m;
    slong guard;
    slong nblocks2;
    ulong * bnd;
    _acb_mzeta_fixacc_struct * acc;   /* sint nres accumulators */
}
series_ctx;

static slong
block_prec(const series_ctx * C, ulong k0)
{
    const side_struct * Sd = C->Sd;
    double V0 = (double) Sd->sh.g * k0 + Sd->sh.r0;
    return taper(Sd->target + Sd->cbits, -Sd->gamd * V0 * 1.4426950408889634) + C->guard;
}

static void
series_worker(slong b, void * vctx)
{
    series_ctx * C = vctx;
    side_struct * Sd = C->Sd;
    ulong m = C->m, N = Sd->sh.N, g = Sd->sh.g, r0 = Sd->sh.r0;
    ulong j0 = C->bnd[b], j1 = C->bnd[b + 1], j, i;
    slong fn_alloc = 0, jj, ri, P = Sd->sint * C->nres, np, k;
    arb_t Q, s, t;
    ulong * facc = NULL, * fwork = NULL, * num, * den, * Vp, * cnt, * errs;
    slong * which;
    int * limbs;
    void * chunks;

    if (j0 >= j1)
        return;

    cnt = flint_malloc(sizeof(ulong) * m * N);

    /* the sums (jj, r), jj = 1, ..., sint, over the same powers of q, are
       formed together (_acb_mzeta_fixed_ratsum_multi reads each table
       entry once for all of them) */
    chunks = flint_malloc(_acb_mzeta_fixed_ratsum_chunks_size(m) * P);
    num = flint_malloc(sizeof(ulong) * m * P);
    den = flint_malloc(sizeof(ulong) * m * P);
    Vp = flint_malloc(sizeof(ulong) * m);
    errs = flint_malloc(sizeof(ulong) * P);
    which = flint_malloc(sizeof(slong) * P);
    limbs = flint_malloc(sizeof(int) * P);
    arb_init(Q);
    arb_init(s);
    arb_init(t);

    /* Q = exp(-gam V(j0 m)) */
    {
        slong p = block_prec(C, j0 * m);
        arb_mul_ui(t, Sd->gam, g * (j0 * m) + r0, p + 16);
        arb_neg(t, t);
        arb_exp(Q, t, p);
    }

    for (j = j0; j < j1; j++)
    {
        ulong k0 = j * m, k1 = FLINT_MIN(Sd->len, k0 + m), len = k1 - k0;
        slong p = block_prec(C, k0);
        slong fn = FLINT_MIN((p + FLINT_BITS - 1) / FLINT_BITS + 1, C->FN);
        slong stride = C->FN + 1, sl;
        const ulong * ftab = C->ftab + (C->FN - fn);
        ulong ftab_err = (fn == C->FN) ? C->ftab_err : 2;

        if (j > j0)
        {
            arb_set_round(t, C->qm, p);
            arb_mul(Q, Q, t, p);
        }

        /* the counts of this block */
        _acb_mzeta_shells_block(cnt, &Sd->sh, k0, len);

        for (i = 0; i < len; i++)
            Vp[i] = 1;

        /* the sums with some nonzero count: numerators and denominators */
        np = 0;
        for (jj = 1; jj <= Sd->sint; jj++)
        {
            for (i = 0; i < len; i++)
            {
                ulong V = g * (k0 + i) + r0;
                if (V != 0)
                    Vp[i] *= V;    /* V^jj < 2^(FLINT_BITS - 2) (checked) */
            }

            for (ri = 0; ri < C->nres; ri++)
            {
                ulong r = C->res[ri], * nu = num + np * len, * de = den + np * len;
                int any = 0;

                for (i = 0; i < len; i++)
                {
                    ulong V = g * (k0 + i) + r0;
                    if (V == 0)
                    {
                        nu[i] = 0;
                        de[i] = 1;
                    }
                    else
                    {
                        nu[i] = cnt[i * N + r];
                        de[i] = Vp[i];
                        any |= (nu[i] != 0);
                    }
                }
                if (!any)
                    continue;

                limbs[np] = _acb_mzeta_fixed_ratsum_choose_limbs(nu, de, len);
                which[np] = (jj - 1) * C->nres + ri;
                np++;
            }
        }

        if (np == 0)
            continue;

        if (fn > fn_alloc)
        {
            fn_alloc = fn;
            facc = flint_realloc(facc, sizeof(ulong) * (fn + 3) * P);
            fwork = flint_realloc(fwork, sizeof(ulong) * _acb_mzeta_fixed_ratsum_multi_work_size(fn, P));
        }

        _acb_mzeta_fixed_ratsum_multi(facc, errs, ftab, stride, fn, num, den, len, np,
            ftab_err, limbs, chunks, fwork);

        for (k = 0; k < np; k++)
        {
            ulong * bs = facc + k * (fn + 3);

            sl = fn + 3;
            while (sl > 0 && bs[sl - 1] == 0)
                sl--;
            if (sl == 0)
                arf_zero(arb_midref(t));
            else
                arf_set_mpn(arb_midref(t), bs, sl, 0);
            arf_mul_2exp_si(arb_midref(t), arb_midref(t), -FLINT_BITS * fn);
            mag_set_ui(arb_radref(t), errs[k]);
            mag_mul_2exp_si(arb_radref(t), arb_radref(t), -FLINT_BITS * fn);

            arb_mul(s, t, Q, p);
            _acb_mzeta_fixacc_add(C->acc + which[k], s, 0);
        }
    }

    arb_clear(Q);
    arb_clear(s);
    arb_clear(t);
    flint_free(facc);
    flint_free(fwork);
    flint_free(chunks);
    flint_free(num);
    flint_free(den);
    flint_free(Vp);
    flint_free(cnt);
    flint_free(errs);
    flint_free(which);
    flint_free(limbs);
}

static void
series_side(side_struct * Sd, slong nthreads)
{
    series_ctx C;
    ulong N = Sd->sh.N, r, nblk;
    slong i, jj, nacc;
    _acb_mzeta_fixacc_struct * acc;
    arb_t q, qm, cur, cj;
    acb_t t;

    arb_init(q);
    arb_init(qm);
    arb_init(cur);
    arb_init(cj);
    acb_init(t);

    C.Sd = Sd;
    C.res = flint_malloc(sizeof(ulong) * N);
    C.nres = 0;
    {
        char * mask = flint_malloc(N);
        _acb_mzeta_shells_residues(mask, &Sd->sh);
        for (r = 0; r < N; r++)
            if (mask[r])
                C.res[C.nres++] = r;
        flint_free(mask);
    }

    nacc = Sd->sint * C.nres;
    acc = flint_malloc(sizeof(_acb_mzeta_fixacc_struct) * FLINT_MAX(nacc, 1));
    for (i = 0; i < nacc; i++)
        _acb_mzeta_fixacc_init(acc + i, Sd->target + Sd->cbits + 32);
    C.acc = acc;

    C.m = FLINT_MAX(1, FLINT_MIN(n_sqrt(2 * Sd->len), 1024));
    nblk = (Sd->len + C.m - 1) / C.m;
    C.guard = 4 + FLINT_BIT_COUNT(nblk) + FLINT_BIT_COUNT(Sd->sint * C.nres);

    /* the shared table of q^i, q = exp(-gam g), i <= m */
    {
        slong pmax = taper(Sd->target + Sd->cbits, 0) + C.guard + 16, wp;
        fmpz_t z;
        mag_t rr;
        ulong * ft;

        fmpz_init(z);
        mag_init(rr);

        C.FN = (pmax + FLINT_BITS - 1) / FLINT_BITS + 1;
        wp = FLINT_BITS * C.FN + 32 + 2 * FLINT_BIT_COUNT(C.m);
        /* gam = pi^2 / (lambda Pi) is recomputed at this precision by the
           caller (Sd->gam has at least wp bits) */
        arb_mul_ui(q, Sd->gam, Sd->sh.g, wp);
        arb_neg(q, q);
        arb_exp(q, q, wp);
        ft = flint_malloc(sizeof(ulong) * (C.m + 1) * (C.FN + 1));
        C.ftab_err = 0;
        arb_one(cur);
        for (i = 0; i <= (slong) C.m; i++)
        {
            arf_get_fmpz_fixed_si(z, arb_midref(cur), -FLINT_BITS * C.FN);
            fmpz_get_ui_array(ft + i * (C.FN + 1), C.FN + 1, z);
            mag_mul_2exp_si(rr, arb_radref(cur), FLINT_BITS * C.FN);
            if (mag_cmp_2exp_si(rr, 40) >= 0)
                flint_throw(FLINT_ERROR, "epstein_shell: q table too inaccurate\n");
            C.ftab_err = FLINT_MAX(C.ftab_err, (ulong) ceil(mag_get_d(rr)) + 1);
            if (i < (slong) C.m)
                arb_mul(cur, cur, q, wp);
        }
        C.ftab = ft;
        arb_set_round(qm, cur, pmax + 16);
        C.qm = qm;

        fmpz_clear(z);
        mag_clear(rr);
    }

    /* tasks of equal estimated cost */
    {
        double * cum = flint_malloc(sizeof(double) * (nblk + 1)), tot;
        ulong jb, bb;

        cum[0] = 0.0;
        for (jb = 0; jb < nblk; jb++)
        {
            ulong k0 = jb * C.m, len = FLINT_MIN(Sd->len, k0 + C.m) - k0;
            cum[jb + 1] = cum[jb] + (double) block_prec(&C, k0) * (len * Sd->sint * C.nres + 8.0);
        }
        tot = cum[nblk];
        C.nblocks2 = FLINT_MAX(1, FLINT_MIN((slong) nblk / 16, 1024));
        C.bnd = flint_malloc(sizeof(ulong) * (C.nblocks2 + 1));
        C.bnd[0] = 0;
        jb = 0;
        for (bb = 1; bb < (ulong) C.nblocks2; bb++)
        {
            double target_cost = tot * bb / C.nblocks2;
            while (jb < nblk && cum[jb] < target_cost)
                jb++;
            C.bnd[bb] = FLINT_MAX(jb, C.bnd[bb - 1]);
        }
        C.bnd[C.nblocks2] = nblk;
        flint_free(cum);
    }

    flint_parallel_do(series_worker, &C, C.nblocks2, nthreads, FLINT_PARALLEL_STRIDED);

    /* sum = sum_r chi_r sum_j (s-1)!/(s-j)! gam^-j S_{j,r} */
    acb_zero(Sd->sum);
    {
        slong wp = Sd->target + Sd->cbits + 32;
        arb_t ig, f, Sj;
        arb_init(ig);
        arb_init(f);
        arb_init(Sj);
        arb_inv(ig, Sd->gam, wp);
        for (i = 0; i < C.nres; i++)
        {
            arb_zero(cj);
            arb_one(f);     /* (s-1)!/(s-j)! gam^-j */
            for (jj = 1; jj <= Sd->sint; jj++)
            {
                arb_mul(f, f, ig, wp);
                if (jj > 1)
                    arb_mul_ui(f, f, Sd->sint - jj + 1, wp);
                _acb_mzeta_fixacc_get(Sj, acc + (jj - 1) * C.nres + i);
                arb_addmul(cj, f, Sj, wp);
            }
            acb_mul_arb(t, Sd->chi + C.res[i], cj, wp);
            acb_add(Sd->sum, Sd->sum, t, wp);
        }
        arb_clear(ig);
        arb_clear(f);
        arb_clear(Sj);
    }

    for (i = 0; i < nacc; i++)
        _acb_mzeta_fixacc_clear(acc + i);
    flint_free(acc);
    flint_free((ulong *) C.ftab);
    flint_free(C.bnd);
    flint_free(C.res);
    arb_clear(q);
    arb_clear(qm);
    arb_clear(cur);
    arb_clear(cj);
    acb_clear(t);
}

/* Setup *********************************************************************/

static void
side_init(side_struct * Sd)
{
    Sd->sh.d = 0;
    fmpq_init(Sd->a);
    arb_init(Sd->gam);
    acb_init(Sd->sum);
    mag_init(Sd->tail);
    Sd->chi = NULL;
}

static void
side_clear(side_struct * Sd)
{
    _acb_mzeta_shells_clear(&Sd->sh);
    fmpq_clear(Sd->a);
    arb_clear(Sd->gam);
    acb_clear(Sd->sum);
    mag_clear(Sd->tail);
    if (Sd->chi != NULL)
        _acb_vec_clear(Sd->chi, Sd->sh.N > 0 ? Sd->sh.N : 1);
}

/* the truncation of a side for |error| <= 2^-bits: len and the rigorous
   tail bound; F_side the form, gam the scale (t = gam S^2 F(n - sigma)) */
static int
side_truncate(side_struct * Sd, slong bits, slong prec)
{
    slong d = Sd->sh.d, i, j;
    arb_mat_t M, Cm;
    arb_t K, t;
    double v, T, Tlo, ad = fmpq_get_d(Sd->a);
    int ok;

    arb_mat_init(M, d, d);
    arb_mat_init(Cm, d, d);
    arb_init(K);
    arb_init(t);

    /* v >= max(0, a - 1), an exact double */
    v = FLINT_MAX(0.0, ceil(ad) - 1.0);

    for (i = 0; i < d; i++)
        for (j = 0; j < d; j++)
        {
            arb_mul_fmpz(arb_mat_entry(M, i, j), Sd->gam, fmpz_mat_entry(Sd->sh.F, i, j), prec);
            arb_mul_ui(arb_mat_entry(M, i, j), arb_mat_entry(M, i, j), Sd->sh.S * Sd->sh.S, prec);
        }
    ok = _acb_mzeta_upper_cho(Cm, M, prec);
    if (ok)
    {
        _acb_mzeta_theta_tail_const(K, Cm, prec);
        T = _acb_mzeta_epstein_choose_T(K, d, v, bits);
        ok = (T < 1e12);
    }
    if (ok)
    {
        double L = (T / Sd->gamd - Sd->sh.r0) / Sd->sh.g + 2.0;
        /* the shell values g len + r0 fit in a ulong */
        ok = (L < 1e12) && ((double) Sd->sh.g * L + Sd->sh.r0 < ldexp(1.0, FLINT_BITS - 3));
        if (ok)
        {
            Sd->len = (ulong) FLINT_MAX(1.0, L);
            /* the first omitted value gam (g len + r0), bounded below */
            arb_mul_ui(t, Sd->gam, Sd->sh.g * Sd->len + Sd->sh.r0, prec);
            Tlo = arf_get_d(arb_midref(t), ARF_RND_DOWN);
            Tlo = Tlo * (1 - 1e-12) - mag_get_d(arb_radref(t));
            _acb_mzeta_epstein_tail_bound(Sd->tail, K, d, Tlo, v);
            ok = mag_is_finite(Sd->tail);
        }
    }

    arb_mat_clear(M);
    arb_mat_clear(Cm);
    arb_clear(K);
    arb_clear(t);
    return ok;
}

static int
fmpq_vec_is_int(const fmpq * x, slong d)
{
    slong i;
    for (i = 0; i < d; i++)
        if (!fmpz_is_one(fmpq_denref(x + i)))
            return 0;
    return 1;
}

static int
is_pos_int(const fmpq_t a)
{
    return fmpz_is_one(fmpq_denref(a)) && fmpz_sgn(fmpq_numref(a)) > 0
        && fmpz_cmp_ui(fmpq_numref(a), 16) <= 0;
}

/* Z_F(u; v) for integral positive definite F (no scaling); returns 0 if
   not applicable */
static int
epstein_shell_int(acb_t res, const fmpq_t nu, const fmpz_mat_t F,
    const fmpq * u, const fmpq * v, slong prec)
{
    slong d = fmpz_mat_nrows(F), i, wp, nthreads = flint_get_num_threads();
    side_struct sides[2], * SD = sides, * SR = sides + 1, * SB, * SO;
    fmpz_mat_t Finv, FR;
    fmpz_t den, h, det;
    fmpq * mv;
    fmpq_t s, b, t;
    arb_t pi, lam, kappa, x, y, PD, PR;
    acb_t z, w, ph;
    double lamd, Pi, cB, cO, digits, Td;
    int ok = 1, qside;

    /* nu = 0, -2, -4, ...: 1/Gamma(s) = 0 and only the singular term
       -e(-v.u) kappa^s / Gamma(s + 1) can survive (s = 0, u in Z^d) */
    if (fmpz_is_one(fmpq_denref(nu)) && fmpz_sgn(fmpq_numref(nu)) <= 0 && fmpz_is_even(fmpq_numref(nu)))
    {
        if (fmpz_is_zero(fmpq_numref(nu)) && fmpq_vec_is_int(u, d))
        {
            fmpq_t c;
            fmpq_init(c);
            for (i = 0; i < d; i++)
                fmpq_addmul(c, u + i, v + i);
            fmpq_mul_si(c, c, -2);
            arb_sin_cos_pi_fmpq(acb_imagref(res), acb_realref(res), c, prec);
            acb_neg(res, res);
            fmpq_clear(c);
        }
        else
            acb_zero(res);
        return 1;
    }

    side_init(SD);
    side_init(SR);
    fmpz_mat_init(Finv, d, d);
    fmpz_mat_init(FR, d, d);
    fmpz_init(den);
    fmpz_init(h);
    fmpz_init(det);
    mv = _fmpq_vec_init(d);
    fmpq_init(s);
    fmpq_init(b);
    fmpq_init(t);
    arb_init(pi);
    arb_init(lam);
    arb_init(kappa);
    arb_init(x);
    arb_init(y);
    arb_init(PD);
    arb_init(PR);
    acb_init(z);
    acb_init(w);
    acb_init(ph);

    /* s = nu/2, b = d/2 - s */
    fmpq_div_2exp(s, nu, 1);
    fmpq_set_si(b, d, 2);
    fmpq_sub(b, b, s);
    if (fmpz_bits(fmpq_denref(s)) > 10 || fmpz_bits(fmpq_numref(s)) > 20)
    {
        ok = 0;
        goto cleanup;
    }

    /* F^-1 = FR / delta */
    fmpz_mat_det(det, F);
    if (fmpz_sgn(det) <= 0 || !fmpz_mat_inv(Finv, den, F))
    {
        ok = 0;
        goto cleanup;
    }
    fmpz_set(h, den);
    for (i = 0; i < d * d; i++)
        fmpz_gcd(h, h, Finv->entries + i);
    if (fmpz_sgn(den) < 0)
        fmpz_neg(h, h);
    fmpz_mat_scalar_divexact_fmpz(FR, Finv, h);
    fmpz_divexact(den, den, h);   /* delta > 0 */

    for (i = 0; i < d; i++)
        fmpq_neg(mv + i, v + i);

    if (!_acb_mzeta_shells_init(&SD->sh, F, u, v) || !_acb_mzeta_shells_init(&SR->sh, FR, mv, u))
    {
        ok = 0;
        goto cleanup;
    }
    fmpq_set(SD->a, s);
    fmpq_set(SR->a, b);
    SD->skip0 = fmpq_vec_is_int(u, d);
    SR->skip0 = fmpq_vec_is_int(v, d);

    /* which side is B */
    if (is_pos_int(b) && (!is_pos_int(s) || SR->sh.N * fmpz_get_ui(fmpq_numref(b)) <= SD->sh.N * fmpz_get_ui(fmpq_numref(s))))
    {
        SB = SD; SO = SR; SO->kind = SIDE_QSERIES;
    }
    else if (is_pos_int(s))
    {
        SB = SR; SO = SD; SO->kind = SIDE_QSERIES;
    }
    else
    {
        SB = SD; SO = SR; SO->kind = SIDE_GENERIC;
    }
    SB->kind = SIDE_BSPLIT;
    qside = (SO->kind == SIDE_QSERIES);
    SB->sint = SO->sint = 0;
    if (qside)
        SO->sint = fmpz_get_si(fmpq_numref(SO->a));

    /* lambda: balance ~T/(lambda g_B) incomplete gamma values of relative
       cost cB against ~T lambda Pi/(pi^2 g_O) terms of cost cO; cB/cO
       for the power series as measured for the Madelung constant */
    digits = prec * 0.30102999566398120;
    Pi = (double) SD->sh.S * SD->sh.S * SR->sh.S * SR->sh.S * fmpz_get_d(den);
    cB = 0.73 * pow(FLINT_MAX(digits, 1.0), 0.8);
    if (qside)
        cO = (double) SO->sint * SO->sh.N;
    else
        cO = cB * 1.5;
    lamd = sqrt(cB / cO * 9.869604401089358 * SO->sh.g / (SB->sh.g * Pi));
    lamd = FLINT_MAX(lamd, 1e-6);
    lamd = FLINT_MIN(lamd, 1e6);

    wp = prec + 32;
    Td = wp * 0.6931471805599453 + 20;

    /* shrink lambda while the power series side would need too many
       counts */
    for (i = 0; i < 64; i++)
    {
        double nO = Td * lamd * Pi / (9.869604401089358 * SO->sh.g) + 2;
        if ((nO < ldexp(1.0, FLINT_BITS - 3) && _acb_mzeta_shells_cost(&SO->sh, (ulong) nO) <= 2e9) || lamd < 1e-6)
            break;
        lamd *= 0.7;
    }

    /* dyadic lambda with 3 bits */
    {
        int e;
        double mfr = frexp(lamd, &e);
        mfr = ceil(mfr * 8) / 8;
        lamd = ldexp(mfr, e);
    }
    arb_set_d(lam, lamd);

    /* precisions: absolute accuracy 2^-target for the side sums, enough
       for the prefactors */
    arb_const_pi(pi, wp + 64);
    if (SB == SD)
    {
        arb_mul_ui(kappa, lam, SD->sh.S * SD->sh.S, ARF_PREC_EXACT);
    }
    else
    {
        arb_sqr(kappa, pi, wp + 64);
        arb_div_fmpz(kappa, kappa, den, wp + 64);
        arb_div_ui(kappa, kappa, SR->sh.S * SR->sh.S, wp + 64);
        arb_div(kappa, kappa, lam, wp + 64);
    }

    /* PD = kappa^s / Gamma(s), PR = pi^(d/2) det^(-1/2) kappa^(s-d/2) / Gamma(s) */
    {
        slong pp = 64;
        arb_t r, k2;
        arb_init(r);
        arb_init(k2);
        arb_set_fmpq(x, s, pp);
        arb_rgamma(r, x, pp);
        arb_pow(PD, kappa, x, pp);
        arb_mul(PD, PD, r, pp);
        fmpq_set_si(t, d, 2);
        fmpq_sub(t, s, t);
        arb_set_fmpq(y, t, pp);
        arb_pow(PR, kappa, y, pp);
        arb_mul(PR, PR, r, pp);
        arb_set_fmpz(k2, det);
        arb_rsqrt(k2, k2, pp);
        arb_mul(PR, PR, k2, pp);
        arb_set_si(k2, d);
        arb_mul_2exp_si(k2, k2, -1);
        arb_pow(k2, pi, k2, pp);
        arb_mul(PR, PR, k2, pp);
        SD->target = wp + FLINT_MAX(0, arf_abs_bound_lt_2exp_si(arb_midref(PD)));
        SR->target = wp + FLINT_MAX(0, arf_abs_bound_lt_2exp_si(arb_midref(PR)));
        if (!arb_is_finite(PD) || !arb_is_finite(PR))
        {
            /* 1/Gamma(s) = 0 handled below; sizes from kappa only */
            SD->target = SR->target = wp + 64;
        }
        arb_clear(r);
        arb_clear(k2);
    }

    /* scales: gam_B = lambda, gam_O = pi^2 / (lambda Pi) */
    {
        slong gp = FLINT_BITS * ((SO->target + 64 + 64 * 4) / FLINT_BITS + 4) + 128;
        arb_t pp;
        arb_init(pp);
        arb_set(SB->gam, lam);
        arb_const_pi(pp, gp);
        arb_sqr(pp, pp, gp);
        arb_div(pp, pp, lam, gp);
        arb_div_ui(pp, pp, SD->sh.S * SD->sh.S, gp);
        arb_div_ui(pp, pp, SR->sh.S * SR->sh.S, gp);
        arb_div_fmpz(pp, pp, den, gp);
        arb_swap(SO->gam, pp);
        arb_clear(pp);
    }
    SB->gamd = arf_get_d(arb_midref(SB->gam), ARF_RND_NEAR);
    SO->gamd = arf_get_d(arb_midref(SO->gam), ARF_RND_NEAR);

    /* truncation */
    for (i = 0; i < 2 && ok; i++)
    {
        side_struct * Sd = sides + i;
        ok = side_truncate(Sd, Sd->target + 4, 64);
        if (ok)
            ok = (_acb_mzeta_shells_cost(&Sd->sh, Sd->len) <= 1e11);
        if (ok && Sd->kind == SIDE_QSERIES)
        {
            /* the denominators V^sint must fit in a limb; else the
               values of G one by one (bit-burst) */
            double Vm = (double) Sd->sh.g * Sd->len + Sd->sh.r0;
            if (!(Sd->sint * log2(Vm + 1) < FLINT_BITS - 2))
                Sd->kind = SIDE_GENERIC;
        }
        /* the shell values g M + r0 (ulong) */
        if (ok && !((double) Sd->sh.g * Sd->len + Sd->sh.r0 < ldexp(1.0, FLINT_BITS - 3)))
            ok = 0;
        if (ok && Sd->kind != SIDE_QSERIES)
            ok = (Sd->len < 1000000);
    }
    if (!ok)
        goto cleanup;

    /* counts and characters */
    for (i = 0; i < 2 && ok; i++)
    {
        side_struct * Sd = sides + i;
        ulong r;
        ok = _acb_mzeta_shells_compute(&Sd->sh, Sd->len);
        if (!ok)
            break;
        Sd->cbits = FLINT_BIT_COUNT(_acb_mzeta_shells_max_count(&Sd->sh));
        Sd->chi = _acb_vec_init(Sd->sh.N);
        for (r = 0; r < Sd->sh.N; r++)
        {
            fmpq_set_si(t, -2 * (slong) r, Sd->sh.N);
            arb_sin_cos_pi_fmpq(acb_imagref(Sd->chi + r), acb_realref(Sd->chi + r), t, Sd->target + 64);
        }
    }
    if (!ok)
        goto cleanup;

    /* the sums */
    gamma_side(SB, lam, nthreads);
    if (SO->kind == SIDE_QSERIES)
        series_side(SO, nthreads);
    else
        gamma_side(SO, lam, nthreads);
    acb_add_error_mag(SD->sum, SD->tail);
    acb_add_error_mag(SR->sum, SR->tail);

    /* assemble */
    {
        slong p = FLINT_MAX(SD->target, SR->target) + 32;
        arb_t r, k2;
        arb_init(r);
        arb_init(k2);

        arb_const_pi(pi, p);
        if (SB == SR)
        {
            arb_sqr(kappa, pi, p);
            arb_div_fmpz(kappa, kappa, den, p);
            arb_div_ui(kappa, kappa, SR->sh.S * SR->sh.S, p);
            arb_div(kappa, kappa, lam, p);
        }
        arb_set_fmpq(x, s, p);
        arb_rgamma(r, x, p);
        arb_pow(PD, kappa, x, p);
        fmpq_set_si(t, d, 2);
        fmpq_sub(t, s, t);
        arb_set_fmpq(y, t, p);
        arb_pow(PR, kappa, y, p);
        arb_set_fmpz(k2, det);
        arb_rsqrt(k2, k2, p);
        arb_mul(PR, PR, k2, p);
        arb_set_si(k2, d);
        arb_mul_2exp_si(k2, k2, -1);
        arb_pow(k2, pi, k2, p);
        arb_mul(PR, PR, k2, p);       /* without 1/Gamma(s) */

        /* e(-u.v) */
        fmpq_zero(t);
        for (i = 0; i < d; i++)
            fmpq_addmul(t, u + i, v + i);
        fmpq_mul_si(t, t, -2);
        arb_sin_cos_pi_fmpq(acb_imagref(ph), acb_realref(ph), t, p);

        acb_mul_arb(z, SD->sum, PD, p);
        acb_mul(w, SR->sum, ph, p);
        acb_mul_arb(w, w, PR, p);
        acb_add(z, z, w, p);

        /* singular term of v in Z^d: PR / (s - d/2) */
        if (SR->skip0)
        {
            fmpq_set_si(t, d, 2);
            fmpq_sub(t, s, t);
            if (fmpq_is_zero(t))
                acb_indeterminate(z);
            else
            {
                arb_set_fmpq(y, t, p);
                arb_div(y, PR, y, p);
                arb_add(acb_realref(z), acb_realref(z), y, p);
            }
        }
        acb_mul_arb(z, z, r, p);

        /* singular term of u in Z^d: -e(-v.u) kappa^s / Gamma(s + 1) (the
           omitted point n = u carries the character e(-v.n));
           ph = e(-u.v) */
        if (SD->skip0)
        {
            arb_add_ui(y, x, 1, p);
            arb_rgamma(y, y, p);
            arb_mul(y, y, PD, p);
            acb_mul_arb(w, ph, y, p);
            acb_sub(z, z, w, p);
        }

        acb_set_round(res, z, prec);

        /* all characters +-1 and e(-u.v) = +-1: every term is real */
        if (SD->sh.N <= 2 && SR->sh.N <= 2)
        {
            fmpq_zero(t);
            for (i = 0; i < d; i++)
                fmpq_addmul(t, u + i, v + i);
            fmpq_mul_2exp(t, t, 1);
            if (fmpz_is_one(fmpq_denref(t)))
                arb_zero(acb_imagref(res));
        }
        arb_clear(r);
        arb_clear(k2);
    }

cleanup:
    side_clear(SD);
    side_clear(SR);
    fmpz_mat_clear(Finv);
    fmpz_mat_clear(FR);
    fmpz_clear(den);
    fmpz_clear(h);
    fmpz_clear(det);
    _fmpq_vec_clear(mv, d);
    fmpq_clear(s);
    fmpq_clear(b);
    fmpq_clear(t);
    arb_clear(pi);
    arb_clear(lam);
    arb_clear(kappa);
    arb_clear(x);
    arb_clear(y);
    arb_clear(PD);
    arb_clear(PR);
    acb_clear(z);
    acb_clear(w);
    acb_clear(ph);

    if (!ok)
        acb_indeterminate(res);
    return ok;
}

/* Public interfaces ***********************************************************/

int
acb_mzeta_epstein_gram_shell(acb_t res, const acb_t nu, const fmpq_mat_t G,
    const fmpq * x, const fmpq * y, int flags, slong prec)
{
    slong d = fmpq_mat_nrows(G), i, j;
    fmpq_t nuq, c;
    fmpz_t L;
    fmpz_mat_t F;
    int ok;

    if (d < 1 || d > 8 || fmpq_mat_ncols(G) != d || (flags & ACB_MZETA_EPSTEIN_REGULARIZED)
        || !arb_is_zero(acb_imagref(nu)) || !arb_is_exact(acb_realref(nu))
        || arf_bits(arb_midref(acb_realref(nu))) > 20
        || arf_cmpabs_2exp_si(arb_midref(acb_realref(nu)), 10) > 0)
    {
        acb_indeterminate(res);
        return 0;
    }

    /* G = (g/L) F with F integral and primitive (content g);
       Z_G = (L/g)^(nu/2) Z_F */
    fmpz_init(L);
    fmpq_init(nuq);
    fmpq_init(c);
    fmpz_mat_init(F, d, d);

    fmpz_one(L);
    /* row by row: the denominators of a row are at stride 2 */
    for (i = 0; i < d; i++)
        _fmpz_vec_lcm_strided(L, fmpq_mat_entry_den(G, i, 0), 2, d, L);
    for (i = 0; i < d; i++)
        for (j = 0; j < d; j++)
        {
            fmpz_divexact(fmpz_mat_entry(F, i, j), L, fmpq_mat_entry_den(G, i, j));
            fmpz_mul(fmpz_mat_entry(F, i, j), fmpz_mat_entry(F, i, j), fmpq_mat_entry_num(G, i, j));
        }

    {
        fmpz_t g;
        fmpz_init(g);
        _fmpz_vec_content(g, F->entries, d * d);
        if (!fmpz_is_zero(g))
        {
            for (i = 0; i < d; i++)
                _fmpz_vec_scalar_divexact_fmpz(fmpz_mat_row(F, i), fmpz_mat_row(F, i), d, g);
        }
        fmpq_set_fmpz_frac(c, L, g);     /* L/g */
        fmpz_clear(g);
    }

    ok = fmpz_bits(fmpq_numref(c)) <= 40 && fmpz_bits(fmpq_denref(c)) <= 40;
    for (i = 0; i < d && ok; i++)
        for (j = 0; j < i && ok; j++)
            ok = fmpz_equal(fmpz_mat_entry(F, i, j), fmpz_mat_entry(F, j, i));
    if (ok)
    {
        arf_get_fmpq(nuq, arb_midref(acb_realref(nu)));
        ok = epstein_shell_int(res, nuq, F, x, y, prec + 8);
    }

    if (!ok)
        acb_indeterminate(res);

    if (ok && !fmpq_is_one(c))
    {
        arb_t t;
        fmpq_t e;
        arb_init(t);
        fmpq_init(e);
        arb_set_fmpq(t, c, prec + 16);
        fmpq_div_2exp(e, nuq, 1);
        arb_pow_fmpq(t, t, e, prec + 16);
        acb_mul_arb(res, res, t, prec);
        arb_clear(t);
        fmpq_clear(e);
    }
    else if (ok)
        acb_set_round(res, res, prec);

    fmpz_clear(L);
    fmpq_clear(nuq);
    fmpq_clear(c);
    fmpz_mat_clear(F);
    return ok;
}

int
acb_mzeta_epstein_shell(acb_t res, const acb_t nu, const arb_mat_t A,
    arb_srcptr x, arb_srcptr y, int flags, acb_mzeta_cache_t cache, slong prec)
{
    slong d = arb_mat_nrows(A), i, j;
    fmpq_mat_t Aq, G, X, Y, AT;
    fmpq * u, * v;
    int ok = 1;

    if (d < 1 || d > 8 || arb_mat_ncols(A) != d)
    {
        acb_indeterminate(res);
        return 0;
    }
    for (i = 0; i < d && ok; i++)
    {
        ok = arb_is_exact(x + i) && arb_is_exact(y + i);
        for (j = 0; j < d && ok; j++)
            ok = arb_is_exact(arb_mat_entry(A, i, j))
                && (arf_is_zero(arb_midref(arb_mat_entry(A, i, j)))
                    || arf_bits(arb_midref(arb_mat_entry(A, i, j))) <= 32);
    }
    if (!ok)
    {
        acb_indeterminate(res);
        return 0;
    }

    fmpq_mat_init(Aq, d, d);
    fmpq_mat_init(AT, d, d);
    fmpq_mat_init(G, d, d);
    fmpq_mat_init(X, d, 1);
    fmpq_mat_init(Y, d, 1);
    u = _fmpq_vec_init(d);
    v = _fmpq_vec_init(d);

    for (i = 0; i < d; i++)
        for (j = 0; j < d; j++)
            arf_get_fmpq(fmpq_mat_entry(Aq, i, j), arb_midref(arb_mat_entry(A, i, j)));
    fmpq_mat_transpose(AT, Aq);
    fmpq_mat_mul(G, AT, Aq);

    /* x = A u, y = A^-T v */
    for (i = 0; i < d; i++)
        arf_get_fmpq(fmpq_mat_entry(X, i, 0), arb_midref(x + i));
    ok = fmpq_mat_solve(Y, Aq, X);
    if (ok)
    {
        for (i = 0; i < d; i++)
            fmpq_set(u + i, fmpq_mat_entry(Y, i, 0));
        for (i = 0; i < d; i++)
            arf_get_fmpq(fmpq_mat_entry(X, i, 0), arb_midref(y + i));
        fmpq_mat_mul(Y, AT, X);
        for (i = 0; i < d; i++)
            fmpq_set(v + i, fmpq_mat_entry(Y, i, 0));
        ok = acb_mzeta_epstein_gram_shell(res, nu, G, u, v, flags, prec);
    }
    else
        acb_indeterminate(res);

    fmpq_mat_clear(Aq);
    fmpq_mat_clear(AT);
    fmpq_mat_clear(G);
    fmpq_mat_clear(X);
    fmpq_mat_clear(Y);
    _fmpq_vec_clear(u, d);
    _fmpq_vec_clear(v, d);
    (void) cache;
    return ok;
}
