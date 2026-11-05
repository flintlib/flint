/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <string.h>
#include "arb_mat.h"
#include "gr.h"
#include "gr_vec.h"
#include "dfloat.h"
#include "double_interval.h"
#include "impl.h"

/* The Crandall sums of acb_mzeta_epstein_crandall in double-precision ball
   arithmetic (dfloat dNb, n = 1 or 2), for real nu and low precision.

   The lattice points are enumerated line by line (coordinate 0 varying)
   in plain double arithmetic with an enlarged radius (Fincke-Pohst; extra
   points are harmless, and every point with t <= T is included), and
   processed in batches through gr vector operations on dNb balls:

     per line:   z_j = n_j + c_j (j >= 1), y_i = sum_{j >= i} C_ij z_j,
                 gam = sum_{i >= 1} y_i^2, beta = C_00 c_0 + sum_{j >= 1} C_0j z_j,
                 and the phase offset delta = sum_{j >= 1} p_j (n_j + pofs_j)
                 + p_0 pofs_0;
     per point:  y_0 = C_00 n_0 + beta, t = y_0^2 + gam, E = exp(-t),
                 F(t) from a table (_acb_mzeta_ftab_eval: plain dN Horner
                 with an a priori bound), G = E F, theta = p_0 n_0 + delta
                 reduced mod 1 exactly, (cos, sin)(2 pi theta), and the
                 dot products of G with cos and sin.

   Points outside the table (t below 2^-12, i.e. very close to the shift)
   are evaluated with _acb_mzeta_gamma_upper_scaled.  The assembly is the
   one of the arb algorithm (_acb_mzeta_epstein_finish). */

#define BATCH 256
#define ELT(v, i) ((double *) (v) + (i) * es)
/* bound for the lattice coordinates (slong) of the enumeration */
#define DFLOAT_COORD_MAX ((FLINT_BITS == 64) ? 0x1p40 : 0x1p28)
#define TMIN_TABLE 0x1p-12

typedef struct
{
    gr_ctx_struct * R;
    slong es;               /* doubles per element: n + 1 */
    slong d;
    gr_ptr C;               /* d x d */
    gr_ptr c, p, pofs;      /* d each */
    gr_ptr twopi;
    int phase;
    slong qe;               /* dyadic phases: p = a / 2^qe exactly (qe >= 0), else -1 */
    slong pa[8];            /* the numerators a_j (mod 2^qe) */
    gr_ptr qcos, qsin;      /* cos, sin of 2 pi r / 2^qe, r < 2^qe */
    slong * pr;             /* per point: r = a.n mod 2^qe */
    di_t Cdi[8][8];         /* enclosures, for the enumeration */
    di_t cdi[8];
    int Cz[8][8], cz[8];    /* exact zeros */
    const _acb_mzeta_ftab_struct * tab;
    const _acb_mzeta_epstein_problem_struct * P;
    int side;
    acb_srcptr a;
    slong prec;

    /* batch */
    slong npts, nlines;
    slong * ln;             /* line coordinates n_1..n_{d-1}: nlines (d - 1) */
    slong * pn0, * pli;     /* per point: n_0, line index */
    gr_ptr lz, ly, lgam, lbeta, ldelta, ltmp;  /* per line (BATCH each; ly: d BATCH) */
    gr_ptr n0, gbeta, ggam, gdelta, y0, t, E, F, G, th, sn, cs, tmp;
    gr_ptr re, im;          /* accumulators */
    gr_ptr bc, bt;          /* broadcast scratch */
    gr_ptr bC00, bp0, b2pi; /* broadcast C_00, p_0, 2 pi (BATCH each) */
    int status;
}
dsum_struct;

static void
set_si_elt(double * x, slong v, slong es)
{
    slong i;
    x[0] = (double) v;
    for (i = 1; i < es; i++)
        x[i] = 0.0;
}

/* vector-scalar operations as vector-vector operations on a broadcast
   copy of the scalar (the vector-vector kernels of the dNb balls are
   much faster than the scalar variants for n >= 2) */
static void
bcast(dsum_struct * S, gr_srcptr c, slong len)
{
    slong i, sz = S->R->sizeof_elem;
    for (i = 0; i < len; i++)
        memcpy(GR_ENTRY(S->bc, i, sz), c, sz);
}

static int
vadd_s(dsum_struct * S, gr_ptr res, gr_srcptr x, slong len, gr_srcptr c)
{
    bcast(S, c, len);
    return _gr_vec_add(res, x, S->bc, len, S->R);
}

static int
vaddmul_s(dsum_struct * S, gr_ptr res, gr_srcptr x, slong len, gr_srcptr c)
{
    int status;
    bcast(S, c, len);
    status = _gr_vec_mul(S->bt, x, S->bc, len, S->R);
    return status | _gr_vec_add(res, res, S->bt, len, S->R);
}

static void
process_batch(dsum_struct * S)
{
    gr_ctx_struct * R = S->R;
    slong es = S->es, d = S->d, L = S->nlines, N = S->npts, i, j, k;
    int status = GR_SUCCESS;

    if (N == 0)
        return;

    /* line quantities */
    status |= _gr_vec_zero(S->lgam, L, R);
    status |= _gr_vec_zero(S->lbeta, L, R);
    for (j = 1; j < d; j++)
    {
        gr_ptr zj = GR_ENTRY(S->lz, (j - 1) * BATCH, R->sizeof_elem);
        for (k = 0; k < L; k++)
            set_si_elt(ELT(zj, k), S->ln[k * (d - 1) + (j - 1)], es);
    }
    if (S->phase)
    {
        status |= _gr_vec_zero(S->ldelta, L, R);
        for (j = 1; j < d; j++)
        {
            gr_ptr zj = GR_ENTRY(S->lz, (j - 1) * BATCH, R->sizeof_elem);
            /* p_j (n_j + pofs_j) */
            status |= vadd_s(S, S->ltmp, zj, L, GR_ENTRY(S->pofs, j, R->sizeof_elem));
            status |= vaddmul_s(S, S->ldelta, S->ltmp, L, GR_ENTRY(S->p, j, R->sizeof_elem));
        }
        /* + p_0 pofs_0 */
        status |= gr_mul(S->tmp, S->p, S->pofs, R);
        status |= vadd_s(S, S->ldelta, S->ldelta, L, S->tmp);
    }
    for (j = 1; j < d; j++)
    {
        gr_ptr zj = GR_ENTRY(S->lz, (j - 1) * BATCH, R->sizeof_elem);
        status |= vadd_s(S, zj, zj, L, GR_ENTRY(S->c, j, R->sizeof_elem));
    }
    for (i = 1; i < d; i++)
    {
        gr_ptr yi = GR_ENTRY(S->ly, (i - 1) * BATCH, R->sizeof_elem);
        status |= _gr_vec_zero(yi, L, R);
        for (j = i; j < d; j++)
        {
            gr_ptr zj = GR_ENTRY(S->lz, (j - 1) * BATCH, R->sizeof_elem);
            status |= vaddmul_s(S, yi, zj, L, GR_ENTRY(S->C, i * d + j, R->sizeof_elem));
        }
        status |= _gr_vec_mul(S->ltmp, yi, yi, L, R);
        status |= _gr_vec_add(S->lgam, S->lgam, S->ltmp, L, R);
    }
    for (j = 1; j < d; j++)
    {
        gr_ptr zj = GR_ENTRY(S->lz, (j - 1) * BATCH, R->sizeof_elem);
        status |= vaddmul_s(S, S->lbeta, zj, L, GR_ENTRY(S->C, j, R->sizeof_elem));
    }
    status |= gr_mul(S->tmp, S->C, S->c, R);
    status |= vadd_s(S, S->lbeta, S->lbeta, L, S->tmp);

    /* point quantities */
    status |= _gr_vec_gather(S->gbeta, S->lbeta, S->pli, N, R);
    status |= _gr_vec_gather(S->ggam, S->lgam, S->pli, N, R);
    for (k = 0; k < N; k++)
        set_si_elt(ELT(S->n0, k), S->pn0[k], es);
    status |= _gr_vec_mul(S->y0, S->n0, S->bC00, N, R);
    status |= _gr_vec_add(S->y0, S->y0, S->gbeta, N, R);
    status |= _gr_vec_mul(S->t, S->y0, S->y0, N, R);
    status |= _gr_vec_add(S->t, S->t, S->ggam, N, R);
    status |= _gr_vec_neg(S->tmp, S->t, N, R);
    status |= _gr_vec_exp(S->E, S->tmp, N, R);
    _acb_mzeta_ftab_eval(S->F, S->t, N, S->tab);
    status |= _gr_vec_mul(S->G, S->E, S->F, N, R);

    /* omitted points, and points outside the table */
    for (k = 0; k < N; k++)
    {
        double * g = ELT(S->G, k);
        slong pt[8];
        int skip;

        pt[0] = S->pn0[k];
        for (j = 1; j < d; j++)
            pt[j] = S->ln[S->pli[k] * (d - 1) + (j - 1)];
        skip = _acb_mzeta_epstein_skip(S->P, S->side, pt);

        if (skip)
        {
            for (j = 0; j < es; j++)
                g[j] = 0.0;
        }
        else if (!(ELT(S->F, k)[es - 1] < INFINITY))
        {
            /* t below the table (close to the shift) */
            arb_t tt;
            acb_t gv;
            _dfloat_ops_t ops = _dfloat_ops(es - 1, 1);
            arb_init(tt);
            acb_init(gv);
            ops->get_arb(tt, ELT(S->t, k));
            _acb_mzeta_gamma_upper_scaled(gv, S->a, tt, NULL, S->prec + 10);
            ops->set_arb(g, acb_realref(gv));
            arb_clear(tt);
            acb_clear(gv);
        }
    }

    if (S->phase && S->qe >= 0)
    {
        /* dyadic phases: e(-p.n) = e(-r / 2^qe), r = a.n mod 2^qe, from a
           table (exact integer arithmetic; the constant e(-p.pofs) is
           applied at the end) */
        ulong mask = (UWORD(1) << S->qe) - 1;
        for (k = 0; k < N; k++)
        {
            ulong r = (ulong) S->pa[0] * (ulong) S->pn0[k];
            for (j = 1; j < d; j++)
                r += (ulong) S->pa[j] * (ulong) S->ln[S->pli[k] * (d - 1) + (j - 1)];
            S->pr[k] = r & mask;
        }
        if (S->qe == 1)
        {
            /* e(-r/2) = (-1)^r: flip the signs (midpoint only) and sum */
            for (k = 0; k < N; k++)
                if (S->pr[k])
                    for (j = 0; j < es - 1; j++)
                        ELT(S->G, k)[j] = -ELT(S->G, k)[j];
            status |= _gr_vec_sum(S->tmp, S->G, N, R);
            status |= gr_add(S->re, S->re, S->tmp, R);
        }
        else
        {
            status |= _gr_vec_gather(S->cs, S->qcos, S->pr, N, R);
            status |= _gr_vec_gather(S->sn, S->qsin, S->pr, N, R);
            status |= _gr_vec_dot(S->re, S->re, 0, S->G, S->cs, N, R);
            status |= _gr_vec_dot(S->im, S->im, 1, S->G, S->sn, N, R);
        }
    }
    else if (S->phase)
    {
        status |= _gr_vec_gather(S->gdelta, S->ldelta, S->pli, N, R);
        status |= _gr_vec_mul(S->th, S->n0, S->bp0, N, R);
        status |= _gr_vec_add(S->th, S->th, S->gdelta, N, R);
        /* exact reduction modulo 1: subtract an integer near the head
           (|head| < 2^52 here), then renormalize the expansion */
        for (k = 0; k < N; k++)
        {
            double * x = ELT(S->th, k);
            if (!(fabs(x[0]) < 0x1p50) || !(x[es - 1] < 1.0))
            {
                /* not reducible exactly: the whole line (sin, cos give
                   [0 +/- 1]) */
                for (j = 0; j < es - 1; j++)
                    x[j] = 0.0;
                x[es - 1] = INFINITY;
                continue;
            }
            x[0] -= nearbyint(x[0]);
            if (es == 3)
            {
                double h, l;
                if (fabs(x[0]) >= fabs(x[1]))
                    DFLOAT_FAST_TWO_SUM(h, l, x[0], x[1]);
                else
                    DFLOAT_TWO_SUM(h, l, x[0], x[1]);
                x[0] = h;
                x[1] = l;
            }
        }
        status |= _gr_vec_mul(S->th, S->th, S->b2pi, N, R);
        status |= _gr_vec_sin_cos(S->sn, S->cs, S->th, N, R);
        /* e(-theta) = cos - i sin */
        status |= _gr_vec_dot(S->re, S->re, 0, S->G, S->cs, N, R);
        status |= _gr_vec_dot(S->im, S->im, 1, S->G, S->sn, N, R);
    }
    else
    {
        status |= _gr_vec_sum(S->tmp, S->G, N, R);
        status |= gr_add(S->re, S->re, S->tmp, R);
    }

    S->status |= status;
    S->npts = 0;
    S->nlines = 0;
}

/* register the line n_0 in [lo, hi] with the other coordinates n_1.. */
static void
add_line(dsum_struct * S, const slong * n, slong lo, slong hi)
{
    slong d = S->d, j, k;

    while (lo <= hi)
    {
        slong cnt;
        if (S->npts == BATCH || S->nlines == BATCH)
            process_batch(S);
        cnt = FLINT_MIN(hi - lo + 1, BATCH - S->npts);
        for (j = 1; j < d; j++)
            S->ln[S->nlines * (d - 1) + (j - 1)] = n[j];
        for (k = 0; k < cnt; k++)
        {
            S->pn0[S->npts + k] = lo + k;
            S->pli[S->npts + k] = S->nlines;
        }
        S->npts += cnt;
        S->nlines++;
        lo += cnt;
    }
}

/* All n with |C (n + c)|^2 <= T, plus possibly a few more (Fincke-Pohst,
   the coordinates from the last to the first), in outward-rounded double
   interval arithmetic (double_interval) on enclosures of C and c, so that
   no point is missed.  Returns the number of points, or -1 if the ranges
   are not finite or unreasonably large. */
static slong
enumerate(dsum_struct * S, double T)
{
    slong d = S->d, n[8] = {0}, lo[8] = {0}, hi[8] = {0}, total = 0, i, j;
    double rem[8];
    int lev;

    /* n_j + c_j; exact zeros are skipped throughout: the outward rounding
       of double_interval turns 0 into tiny subnormal intervals, which are
       very slow to compute with */
#define ZN(jj) (S->cz[jj] ? di_interval((double) n[jj], (double) n[jj]) : di_fast_add_d(S->cdi[jj], (double) n[jj]))

    /* range of coordinate i given n_{i+1..d-1}: |C_ii n_i + beta| <= sqrt(rem_i),
       beta = C_ii c_i + sum_{j > i} C_ij (n_j + c_j) */
#define SET_RANGE(ii) do { \
        di_t bta = S->cz[ii] ? di_interval(0.0, 0.0) : di_fast_mul(S->Cdi[ii][ii], S->cdi[ii]), r, x, y; \
        for (j = (ii) + 1; j < d; j++) \
            if (!S->Cz[ii][j]) \
                bta = di_fast_add(bta, di_fast_mul(S->Cdi[ii][j], ZN(j))); \
        r = di_interval(0.0, _di_above(sqrt(FLINT_MAX(rem[ii], 0.0)))); \
        x = di_fast_div(di_fast_sub(di_neg(r), bta), S->Cdi[ii][ii]); \
        y = di_fast_div(di_fast_sub(r, bta), S->Cdi[ii][ii]); \
        if (!(x.a > -DFLOAT_COORD_MAX && y.b < DFLOAT_COORD_MAX)) \
            return -1; \
        lo[ii] = (slong) ceil(x.a); \
        hi[ii] = (slong) floor(y.b); \
    } while (0)

    for (i = 0; i < d; i++)
        if (!(S->Cdi[i][i].a > 0.0))
            return -1;

    lev = d - 1;
    rem[lev] = T;
    SET_RANGE(lev);
    n[lev] = lo[lev];
    for (;;)
    {
        if (n[lev] > hi[lev])
        {
            if (lev == d - 1)
                break;
            lev++;
            n[lev]++;
            continue;
        }
        if (lev > 0)
        {
            /* rem_{lev-1} = rem_lev - y^2 (rounded up),
               y = C_ll (n_l + c_l) + sum_{j > l} C_lj (n_j + c_j) */
            di_t yy = di_fast_mul(S->Cdi[lev][lev], ZN(lev));
            for (j = lev + 1; j < d; j++)
                if (!S->Cz[lev][j])
                    yy = di_fast_add(yy, di_fast_mul(S->Cdi[lev][j], ZN(j)));
            rem[lev - 1] = _di_above(rem[lev] - di_fast_sqr(yy).a);
            if (!(rem[lev - 1] < 1e300))
                return -1;
            lev--;
            SET_RANGE(lev);
            n[lev] = lo[lev];
            continue;
        }
        if (lo[0] <= hi[0])
        {
            add_line(S, n, lo[0], hi[0]);
            total += hi[0] - lo[0] + 1;
        }
        if (d == 1)
            break;
        lev = 1;
        n[1]++;
    }
#undef SET_RANGE
#undef ZN

    process_batch(S);
    return total;
}

static int
dfloat_side(acb_t res, const _acb_mzeta_epstein_problem_struct * P, int side,
    const _acb_mzeta_lattice_struct * L, int nn, acb_mzeta_cache_t cache, slong prec)
{
    gr_ctx_t R;
    dsum_struct S;
    slong d = L->d, i, j, es = nn + 1;
    const arb_mat_struct * C = side ? L->Cr : L->Cd;
    const arb_struct * K = side ? L->Kr : L->Kd;
    double v = side ? P->v1 : P->v0, T, a;
    mag_t tail;
    arb_t tmp;
    _dfloat_ops_t ops = _dfloat_ops(nn, 1);
    int ok = 1;

    if (gr_ctx_init_dfloat(R, nn, DFLOAT_BALL | DFLOAT_FAST) != GR_SUCCESS)
        return 0;

    mag_init(tail);
    arb_init(tmp);

    /* the table is for this exact index: nu/2 and (d - nu)/2 need not
       both be doubles when nu is */
    if (arf_bits(arb_midref(acb_realref(side ? P->b : P->a))) > 53)
    {
        mag_clear(tail);
        arb_clear(tmp);
        gr_ctx_clear(R);
        return 0;
    }
    a = arf_get_d(arb_midref(acb_realref(side ? P->b : P->a)), ARF_RND_NEAR);
    T = _acb_mzeta_epstein_choose_T(K, d, v, prec + 8);
    /* (NaN fails) a huge truncation point: the lattice is too skewed for
       this algorithm */
    if (!(T < 1e6))
    {
        mag_clear(tail);
        arb_clear(tmp);
        gr_ctx_clear(R);
        return 0;
    }
    _acb_mzeta_epstein_tail_bound(tail, K, d, T, v);

    S.R = R;
    S.es = es;
    S.d = d;
    S.P = P;
    S.side = side;
    S.a = side ? P->b : P->a;
    S.prec = prec;
    S.status = GR_SUCCESS;
    S.npts = S.nlines = 0;
    S.tab = _acb_mzeta_ftab_get(cache, a, nn, TMIN_TABLE, T + 2.0);

    GR_TMP_INIT_VEC(S.C, d * d, R);
    GR_TMP_INIT_VEC(S.c, d, R);
    GR_TMP_INIT_VEC(S.p, d, R);
    GR_TMP_INIT_VEC(S.pofs, d, R);
    GR_TMP_INIT_VEC(S.twopi, 1, R);
    GR_TMP_INIT_VEC(S.re, 1, R);
    GR_TMP_INIT_VEC(S.im, 1, R);
    GR_TMP_INIT_VEC(S.tmp, BATCH, R);
    GR_TMP_INIT_VEC(S.lz, FLINT_MAX(d - 1, 1) * BATCH, R);
    GR_TMP_INIT_VEC(S.ly, FLINT_MAX(d - 1, 1) * BATCH, R);
    GR_TMP_INIT_VEC(S.lgam, BATCH, R);
    GR_TMP_INIT_VEC(S.lbeta, BATCH, R);
    GR_TMP_INIT_VEC(S.ldelta, BATCH, R);
    GR_TMP_INIT_VEC(S.ltmp, BATCH, R);
    GR_TMP_INIT_VEC(S.n0, BATCH, R);
    GR_TMP_INIT_VEC(S.gbeta, BATCH, R);
    GR_TMP_INIT_VEC(S.ggam, BATCH, R);
    GR_TMP_INIT_VEC(S.gdelta, BATCH, R);
    GR_TMP_INIT_VEC(S.y0, BATCH, R);
    GR_TMP_INIT_VEC(S.t, BATCH, R);
    GR_TMP_INIT_VEC(S.E, BATCH, R);
    GR_TMP_INIT_VEC(S.F, BATCH, R);
    GR_TMP_INIT_VEC(S.G, BATCH, R);
    GR_TMP_INIT_VEC(S.th, BATCH, R);
    GR_TMP_INIT_VEC(S.sn, BATCH, R);
    GR_TMP_INIT_VEC(S.cs, BATCH, R);
    GR_TMP_INIT_VEC(S.bc, BATCH, R);
    GR_TMP_INIT_VEC(S.bC00, BATCH, R);
    GR_TMP_INIT_VEC(S.bp0, BATCH, R);
    GR_TMP_INIT_VEC(S.b2pi, BATCH, R);
    GR_TMP_INIT_VEC(S.bt, BATCH, R);
    S.ln = flint_malloc(BATCH * FLINT_MAX(d - 1, 1) * sizeof(slong));
    S.pn0 = flint_malloc(BATCH * sizeof(slong));
    S.pli = flint_malloc(BATCH * sizeof(slong));

    /* the quadratic form and the phase:
       direct: center c = -u0, phase e(-w0.n) (pofs = 0);
       reciprocal: center c = w0, phase e(-u0.(m + w0)) */
    for (i = 0; i < d; i++)
    {
        for (j = 0; j < d; j++)
        {
            ops->set_arb(ELT(S.C, i * d + j), arb_mat_entry(C, i, j));
            S.Cdi[i][j] = arb_get_di(arb_mat_entry(C, i, j));
            S.Cz[i][j] = arb_is_zero(arb_mat_entry(C, i, j));
        }
        if (side == 0)
        {
            arb_neg(tmp, P->u0 + i);
            ops->set_arb(ELT(S.c, i), tmp);
            ops->set_arb(ELT(S.p, i), P->w0 + i);
            ops->zero(ELT(S.pofs, i));
        }
        else
        {
            ops->set_arb(ELT(S.c, i), P->w0 + i);
            ops->set_arb(ELT(S.p, i), P->u0 + i);
            ops->set_arb(ELT(S.pofs, i), P->w0 + i);
        }
        if (side == 0)
        {
            arb_neg(tmp, P->u0 + i);
            S.cdi[i] = arb_get_di(tmp);
        }
        else
            S.cdi[i] = arb_get_di(P->w0 + i);
        S.cz[i] = arb_is_zero(side ? P->w0 + i : P->u0 + i);
        ok = ok && S.Cdi[i][i].a > 0.0 && S.cdi[i].a > -1.0 && S.cdi[i].b < 1.0;
    }
    /* e(-p.(n + pofs)) = 1 when p = 0 */
    S.phase = !_arb_vec_is_zero(side ? P->u0 : P->w0, d);

    arb_const_pi(tmp, 64 * nn + 64);
    arb_mul_2exp_si(tmp, tmp, 1);
    ops->set_arb(S.twopi, tmp);
    /* dyadic phases p = a / 2^qe, qe <= 6 (e.g. alternating sums) */
    S.qe = -1;
    S.qcos = S.qsin = NULL;
    S.pr = NULL;
    if (S.phase)
    {
        arb_srcptr pv = side ? P->u0 : P->w0;
        slong e, qe = 0;
        int dy = 1;
        for (i = 0; i < d && dy; i++)
        {
            if (!arb_is_exact(pv + i))
                dy = 0;
            else if (!arf_is_zero(arb_midref(pv + i)))
            {
                e = -(slong) fmpz_get_si(ARF_EXPREF(arb_midref(pv + i))) + (slong) arf_bits(arb_midref(pv + i));
                if (!COEFF_IS_MPZ(*ARF_EXPREF(arb_midref(pv + i))) && e <= 6)
                    qe = FLINT_MAX(qe, e);
                else
                    dy = 0;
            }
        }
        if (dy && qe >= 1)
        {
            slong q = WORD(1) << qe, r;
            S.qe = qe;
            for (i = 0; i < d; i++)
            {
                arb_mul_2exp_si(tmp, pv + i, qe);
                S.pa[i] = arf_get_si(arb_midref(tmp), ARF_RND_DOWN) & (q - 1);
            }
            GR_TMP_INIT_VEC(S.qcos, q, R);
            GR_TMP_INIT_VEC(S.qsin, q, R);
            S.pr = flint_malloc(BATCH * sizeof(slong));
            /* the table by the vector kernels: angles 2 pi r 2^-qe */
            for (r = 0; r < q; r++)
            {
                set_si_elt(ELT(S.qsin, r), r, es);
                ops->mul_2exp_si(ELT(S.qsin, r), ELT(S.qsin, r), -qe);
            }
            S.status |= _gr_vec_mul_scalar(S.qsin, S.qsin, q, S.twopi, R);
            S.status |= _gr_vec_sin_cos(S.qsin, S.qcos, S.qsin, q, R);
        }
    }
    bcast(&S, S.C, BATCH);
    memcpy(S.bC00, S.bc, BATCH * R->sizeof_elem);
    bcast(&S, S.p, BATCH);
    memcpy(S.bp0, S.bc, BATCH * R->sizeof_elem);
    bcast(&S, S.twopi, BATCH);
    memcpy(S.b2pi, S.bc, BATCH * R->sizeof_elem);
    ops->zero(S.re);
    ops->zero(S.im);

    if (ok)
    {
        ok = (enumerate(&S, T) >= 0);
        ops->get_arb(acb_realref(res), S.re);
        ops->get_arb(acb_imagref(res), S.im);
        if (S.qe >= 0 && side == 1)
        {
            /* the constant factor e(-p.pofs) = e(-u0.w0) */
            arb_t ang;
            acb_t e;
            arb_init(ang);
            acb_init(e);
            for (i = 0; i < d; i++)
                arb_addmul(ang, P->u0 + i, P->w0 + i, 64 * nn + 64);
            arb_mul_2exp_si(ang, ang, 1);
            arb_neg(ang, ang);
            arb_sin_cos_pi(acb_imagref(e), acb_realref(e), ang, 64 * nn + 64);
            acb_mul(res, res, e, 64 * nn + 64);
            arb_clear(ang);
            acb_clear(e);
        }
        acb_add_error_mag(res, tail);
        ok = ok && (S.status == GR_SUCCESS) && acb_is_finite(res);
    }

    if (S.qe >= 0)
    {
        GR_TMP_CLEAR_VEC(S.qcos, WORD(1) << S.qe, R);
        GR_TMP_CLEAR_VEC(S.qsin, WORD(1) << S.qe, R);
        flint_free(S.pr);
    }
    GR_TMP_CLEAR_VEC(S.C, d * d, R);
    GR_TMP_CLEAR_VEC(S.c, d, R);
    GR_TMP_CLEAR_VEC(S.p, d, R);
    GR_TMP_CLEAR_VEC(S.pofs, d, R);
    GR_TMP_CLEAR_VEC(S.twopi, 1, R);
    GR_TMP_CLEAR_VEC(S.re, 1, R);
    GR_TMP_CLEAR_VEC(S.im, 1, R);
    GR_TMP_CLEAR_VEC(S.tmp, BATCH, R);
    GR_TMP_CLEAR_VEC(S.lz, FLINT_MAX(d - 1, 1) * BATCH, R);
    GR_TMP_CLEAR_VEC(S.ly, FLINT_MAX(d - 1, 1) * BATCH, R);
    GR_TMP_CLEAR_VEC(S.lgam, BATCH, R);
    GR_TMP_CLEAR_VEC(S.lbeta, BATCH, R);
    GR_TMP_CLEAR_VEC(S.ldelta, BATCH, R);
    GR_TMP_CLEAR_VEC(S.ltmp, BATCH, R);
    GR_TMP_CLEAR_VEC(S.n0, BATCH, R);
    GR_TMP_CLEAR_VEC(S.gbeta, BATCH, R);
    GR_TMP_CLEAR_VEC(S.ggam, BATCH, R);
    GR_TMP_CLEAR_VEC(S.gdelta, BATCH, R);
    GR_TMP_CLEAR_VEC(S.y0, BATCH, R);
    GR_TMP_CLEAR_VEC(S.t, BATCH, R);
    GR_TMP_CLEAR_VEC(S.E, BATCH, R);
    GR_TMP_CLEAR_VEC(S.F, BATCH, R);
    GR_TMP_CLEAR_VEC(S.G, BATCH, R);
    GR_TMP_CLEAR_VEC(S.th, BATCH, R);
    GR_TMP_CLEAR_VEC(S.sn, BATCH, R);
    GR_TMP_CLEAR_VEC(S.cs, BATCH, R);
    GR_TMP_CLEAR_VEC(S.bc, BATCH, R);
    GR_TMP_CLEAR_VEC(S.bC00, BATCH, R);
    GR_TMP_CLEAR_VEC(S.bp0, BATCH, R);
    GR_TMP_CLEAR_VEC(S.b2pi, BATCH, R);
    GR_TMP_CLEAR_VEC(S.bt, BATCH, R);
    flint_free(S.ln);
    flint_free(S.pn0);
    flint_free(S.pli);
    mag_clear(tail);
    arb_clear(tmp);
    gr_ctx_clear(R);
    return ok;
}

int
acb_mzeta_epstein_crandall_dfloat(acb_t res, const acb_t nu, const arb_mat_t A,
    arb_srcptr x, arb_srcptr y, int flags, acb_mzeta_cache_t cache, slong prec)
{
    slong d = arb_mat_nrows(A);
    int nn, ok;
    const _acb_mzeta_lattice_struct * L;
    _acb_mzeta_epstein_problem_struct P;
    acb_mzeta_cache_t tmp_cache;
    acb_t Sd, Sr;
    double nud;

    /* applicability: real exact nu of moderate size, d <= 8, low precision */
    if (prec > 95 || d < 1 || d > 8 || arb_mat_ncols(A) != d
        || !arb_is_zero(acb_imagref(nu)) || !arb_is_exact(acb_realref(nu))
        || arf_cmpabs_2exp_si(arb_midref(acb_realref(nu)), 6) > 0
        || !dfloat_is_supported())
    {
        acb_indeterminate(res);
        return 0;
    }
    nud = arf_get_d(arb_midref(acb_realref(nu)), ARF_RND_NEAR);
    {
        /* the table index must be exact as a double */
        arf_t t;
        arf_init(t);
        arf_set_d(t, nud);
        ok = arf_equal(t, arb_midref(acb_realref(nu)));
        arf_clear(t);
        if (!ok)
        {
            acb_indeterminate(res);
            return 0;
        }
    }

    if (cache == NULL)
    {
        acb_mzeta_cache_init(tmp_cache);
        ok = acb_mzeta_epstein_crandall_dfloat(res, nu, A, x, y, flags, tmp_cache, prec);
        acb_mzeta_cache_clear(tmp_cache);
        return ok;
    }

    nn = (prec <= 45) ? 1 : 2;

    L = _acb_mzeta_lattice_get(cache, A, 128);
    if (L == NULL)
    {
        acb_indeterminate(res);
        return 0;
    }

    if (!_acb_mzeta_epstein_problem_init(&P, L, nu, x, y, flags, 192))
    {
        _acb_mzeta_epstein_problem_clear(&P);
        acb_indeterminate(res);
        return 1;
    }

    acb_init(Sd);
    acb_init(Sr);

    ok = dfloat_side(Sd, &P, 0, L, nn, cache, prec)
      && dfloat_side(Sr, &P, 1, L, nn, cache, prec);

    if (ok)
        _acb_mzeta_epstein_finish(res, Sd, Sr, L, &P, prec, 64 * nn + 32);
    else
        acb_indeterminate(res);

    acb_clear(Sd);
    acb_clear(Sr);
    _acb_mzeta_epstein_problem_clear(&P);
    return ok;
}
