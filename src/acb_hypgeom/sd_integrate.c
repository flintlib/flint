/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.2

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Rigorous Gauss-Legendre quadrature of exp(phi(w)) along the
    steepest-descent contours of sd_plan.c.

    The integrand is evaluated as exp(phi(w) - c) with c = phi(w_s) at
    the saddle of each piece of the contour (an exact complex double), so
    that the dominant value is factored out; the pieces are then
    combined, and the factor exp(c) restored, in acb arithmetic.  The
    vertices of the polygons are doubles.  On a chord [a, b] between two
    vertices, w = m + h t, the n-point Gauss-Legendre rule has error at
    most

        (64/15) M |h| / ((rho - 1) rho^(2n-1)),

    M = sup |exp(phi(w) - c)| over the Bernstein ellipse with parameter
    rho = X + sqrt(X^2 - 1); the ellipse is contained in the disk
    |w - m| <= R = |h| X, over which

        Re(phi(m+u) - phi(m)) <= |phi'(m)| R + |phi''(m)| (cosh R - 1)
                                 + |phi'''(m)| (sinh R - R),

    plus r_alpha e^(x+R) + r_beta e^(-x+R) + r_nu (|m| + R) for the input
    balls.  The number of nodes and X are chosen per chord for the
    target error, and the panels of the planned polygon are merged into
    longer chords while this does not increase the estimated number of
    nodes.  Panels on which the integrand is below the tolerance are
    bounded (by the majorants of the planner) instead of integrated, as
    are the tails beyond the polygons.  The input balls are used directly
    in the ball evaluation, so the result encloses the integral on the
    whole input balls; the tolerance is relaxed to what the input balls
    permit.

    A finite starting point (a ball) is joined to its double
    approximation, the first vertex of the contour, by one more chord
    with a ball endpoint; this chord is short, and a rule with few nodes
    is certified on a disk much larger than the chord.

    There are two implementations of the rule: in vectorized dfloat ball
    arithmetic with N <= 4 doubles (the precision chosen so that the
    cancellation in phi(w) - c, of size |alpha| e^x + |beta| e^(-x) +
    |nu| |w|, leaves prec bits), and in acb arithmetic at any precision.
    In both, the precision is lowered where the integrand is small: a
    node where |exp(phi - c)| is 2^(-k) needs k bits less (per chord
    with dfloat, per node with acb).
*/

#include "thread_support.h"
#include "arb_hypgeom.h"
#include "acb.h"
#include "acb_calc/impl.h"
#include "acb_hypgeom/sd.h"
#include "dfloat.h"

/* The rules used are those of the cache of acb_calc_integrate, with
   numbers of nodes in steps of about 2^(1/2), up to 4096 nodes (steps
   of 2^(1/4) would save 8-10% of the nodes, at the price of about
   twice as many rules to precompute). */
#define SDQ_STEPS 23
#define sdq_steps _acb_calc_gl_steps
/* the largest rule with dfloat: 256 nodes */
#define SDQ_DFLOAT_STEPS 15
#define SDQ_DFLOAT_NMAX 256
/* bits of guard beyond prec and the cancellation */
#define SDQ_GUARD 4.0

/* per-N kernels */
#define N 1
#define DB(name) d1b_##name
#define DBV(name) _d1b_##name
#define QUAD(name) sdq_##name##_d1
#include "sd_quad.inc"
#undef N
#undef DB
#undef DBV
#undef QUAD

#define N 2
#define DB(name) d2b_##name
#define DBV(name) _d2b_##name
#define QUAD(name) sdq_##name##_d2
#include "sd_quad.inc"
#undef N
#undef DB
#undef DBV
#undef QUAD

#define N 3
#define DB(name) d3b_##name
#define DBV(name) _d3b_##name
#define QUAD(name) sdq_##name##_d3
#include "sd_quad.inc"
#undef N
#undef DB
#undef DBV
#undef QUAD

#define N 4
#define DB(name) d4b_##name
#define DBV(name) _d4b_##name
#define QUAD(name) sdq_##name##_d4
#include "sd_quad.inc"
#undef N
#undef DB
#undef DBV
#undef QUAD

/* The nodes and weights as arb balls come from the cache of
   acb_calc_integrate (the (n + 1) / 2 nonnegative nodes in decreasing
   order, and their weights, possibly with more precision than asked
   for); the cache here holds all the nodes and weights of the rules
   used with dfloat, converted to arrays of dNb balls. */
typedef struct
{
    void * dx[4][SDQ_DFLOAT_STEPS];
    void * dw[4][SDQ_DFLOAT_STEPS];
}
sdq_cache_struct;

FLINT_TLS_PREFIX sdq_cache_struct * sdq_cache = NULL;

static void
sdq_cleanup(void)
{
    int N, i;

    if (sdq_cache == NULL)
        return;

    for (N = 0; N < 4; N++)
        for (i = 0; i < SDQ_DFLOAT_STEPS; i++)
        {
            flint_free(sdq_cache->dx[N][i]);
            flint_free(sdq_cache->dw[N][i]);
        }
    flint_free(sdq_cache);
    sdq_cache = NULL;
}

static void
sdq_get_nodes_dfloat(const void ** x, const void ** w, int N, int i)
{
    arb_srcptr xa, wa;

    if (sdq_cache == NULL)
    {
        sdq_cache = flint_calloc(1, sizeof(sdq_cache_struct));
        flint_register_cleanup_function(sdq_cleanup);
    }

    if (sdq_cache->dx[N - 1][i] == NULL)
    {
        slong n = sdq_steps[i];
        void * xx, * ww;

        _acb_calc_gl_nodes(&xa, &wa, i, 53 * N + 30);

        xx = flint_malloc(n * (N + 1) * sizeof(double));
        ww = flint_malloc(n * (N + 1) * sizeof(double));
        if (N == 1) sdq_nodes_d1(xx, ww, xa, wa, n);
        else if (N == 2) sdq_nodes_d2(xx, ww, xa, wa, n);
        else if (N == 3) sdq_nodes_d3(xx, ww, xa, wa, n);
        else sdq_nodes_d4(xx, ww, xa, wa, n);
        sdq_cache->dx[N - 1][i] = xx;
        sdq_cache->dw[N - 1][i] = ww;
    }

    *x = sdq_cache->dx[N - 1][i];
    *w = sdq_cache->dw[N - 1][i];
}

typedef union
{
    sdq_params_d1 d1;
    sdq_params_d2 d2;
    sdq_params_d3 d3;
    sdq_params_d4 d4;
}
sdq_params;

typedef struct
{
    const sd_problem * P;
    acb_srcptr alpha, nu, start;
    int sym;
    int N;                   /* the largest dfloat precision, or 0 for acb */
    double bits;             /* prec and guard bits, plus the cancellation */
    slong wp;                /* working precision for acb arithmetic */
    int imax;                /* the largest rule */
    slong budget;            /* of nodes */
    sdq_params par[4];
    int par_set[4];
    void * scratch;
}
sdq_ctx;

/* The smallest rule (its index, or -1) with a certified error below
   exp(logtol) on any chord with half-length at most ah whose
   Bernstein ellipses with parameter X are contained in the disks
   |w - m| <= ah X + r0, and the log of that error bound: the supremum
   of |exp(phi(w) - c)| (Re c = top) over the disk and the input balls
   is bounded as in the header comment, with a rounding margin.  Also
   sets *nreq to the number of nodes required if it could be any real
   number. */
static int
sdq_choose(double * logerr, double * nreq, const sdq_ctx * ctx, cd m, double ah, double r0, double top, double logtol)
{
    static const double Xs[] = { 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0,
        1e2, 1e3, 1e4, 1e6, 1e8, 1e10, 1e12 };
    const sd_problem * P = ctx->P;
    cd phi, d1, d2;
    double am = cd_abs(m), A, B, D, base, bestn = INFINITY, bestneed = 0.0, bestlrho = 1.0;
    int ix, i;

    _sd_phi(&phi, &d1, &d2, P, m);
    A = cd_abs(d2);
    B = cd_abs(cd_add(d1, P->nu));
    D = cd_abs(d1);
    base = phi.re - top + log(64.0 / 15.0) + log(ah) - logtol
        + P->rnu * am + 64.0 * SD_EPS * (P->anu * am + fabs(top) + 1.0);

    for (ix = 0; ix < (int) (sizeof(Xs) / sizeof(double)); ix++)
    {
        double X = Xs[ix], rho = X + sqrt(X * X - 1.0), lrho = log(rho), R = ah * X + r0;
        double sh, chm1, shmr, need, n, size, ball;

        if (R > 300.0)
            break;

        /* cosh R - 1 and sinh R - R, rounded up */
        sh = sinh(0.5 * R);
        chm1 = 2.0 * sh * sh;
        shmr = (R < 0.5) ? R * R * R * 0.17 : sinh(R) - R;
        size = P->aal * exp(m.re + R) + ((P->abe != 0.0) ? P->abe * exp(-m.re + R) : 0.0);
        ball = P->ral * exp(m.re + R) + ((P->abe != 0.0) ? P->rbe * exp(-m.re + R) : 0.0);

        /* (2n-1) log rho >= log(64/15) + log M + log(ah) - log(rho-1) - logtol */
        need = base + (D * R + A * chm1 + B * shmr) * (1.0 + 1e-12)
            + ball + P->rnu * R + 64.0 * SD_EPS * (size + P->anu * R) - log(rho - 1.0);
        n = 0.5 * (need / lrho + 1.0);

        /* (the requirement first decreases and then increases with X) */
        if (!(n < bestn))
            break;
        bestn = n;
        bestneed = need;
        bestlrho = lrho;
    }

    *nreq = fmax(bestn, 1.0);

    for (i = 0; i <= ctx->imax; i++)
        if (sdq_steps[i] >= bestn)
            break;
    if (i > ctx->imax)
        return -1;

    *logerr = bestneed - (2 * sdq_steps[i] - 1) * bestlrho + logtol;
    return i;
}

static int
sdq_choose_chord(double * logerr, double * nreq, const sdq_ctx * ctx, cd a, cd b, double top, double logtol)
{
    return sdq_choose(logerr, nreq, ctx, cd_scale(cd_add(a, b), 0.5),
        0.5 * cd_abs(cd_sub(b, a)), 0.0, top, logtol);
}

/* distance from p to the segment [a, b] */
static double
sdq_segdist(cd p, cd a, cd b)
{
    cd d = cd_sub(b, a), q = cd_sub(p, a);
    double l2 = d.re * d.re + d.im * d.im, t;

    if (l2 == 0.0)
        return cd_abs(q);
    t = (q.re * d.re + q.im * d.im) / l2;
    t = fmin(fmax(t, 0.0), 1.0);
    return cd_abs(cd_sub(q, cd_scale(d, t)));
}

static void
acb_set_cd(acb_t res, cd z)
{
    arb_set_d(acb_realref(res), z.re);
    arb_set_d(acb_imagref(res), z.im);
}

/* The rule sdq_steps[gi] on the chord [a, b] in acb arithmetic.  The
   endpoints are balls; ad and bd are their double approximations, used
   to estimate the size of the integrand at each node, which is then
   evaluated with the corresponding number of bits less. */
static void
sdq_chord_acb(acb_t res, const sdq_ctx * ctx, cd c, const acb_t a, const acb_t b,
    cd ad, cd bd, int gi)
{
    cd md = cd_scale(cd_add(ad, bd), 0.5), hd = cd_scale(cd_sub(bd, ad), 0.5), phi, d1, d2;
    slong n = sdq_steps[gi], k, kk, wp = ctx->wp, p;
    arb_srcptr x, wts;
    acb_t m, h, w, e, t, sum, cc;
    arb_t sh, ch, sn, cs;
    double xd;

    acb_init(m); acb_init(h); acb_init(w); acb_init(e);
    acb_init(t); acb_init(sum); acb_init(cc);
    arb_init(sh); arb_init(ch); arb_init(sn); arb_init(cs);

    _acb_calc_gl_nodes(&x, &wts, gi, wp);

    acb_add(m, a, b, wp);
    acb_mul_2exp_si(m, m, -1);
    acb_sub(h, b, a, wp);
    acb_mul_2exp_si(h, h, -1);
    acb_set_cd(cc, c);

    for (k = 0; k < n; k++)
    {
        kk = (2 * k < n) ? k : n - 1 - k;
        xd = arf_get_d(arb_midref(x + kk), ARF_RND_NEAR);
        if (2 * k >= n)
            xd = -xd;

        _sd_phi(&phi, &d1, &d2, ctx->P, cd_add(md, cd_scale(hd, xd)));
        p = wp + (slong) fmin((phi.re - c.re) / SD_LN2 + 2.0, 0.0);
        p = FLINT_MAX(p, 32);

        /* w = m +- h x (the cached nodes can have a higher precision) */
        arb_set_round(sh, x + kk, p);
        acb_mul_arb(w, h, sh, p);
        if (2 * k >= n)
            acb_sub(w, m, w, p);
        else
            acb_add(w, m, w, p);

        /* E = e^w - e^-w = 2 (sinh x cos y + i cosh x sin y), or e^w */
        arb_sin_cos(sn, cs, acb_imagref(w), p);
        if (ctx->sym)
        {
            arb_sinh_cosh(sh, ch, acb_realref(w), p);
            arb_mul(acb_realref(e), sh, cs, p);
            arb_mul(acb_imagref(e), ch, sn, p);
            acb_mul_2exp_si(e, e, 1);
        }
        else
        {
            arb_exp(ch, acb_realref(w), p);
            arb_mul(acb_realref(e), ch, cs, p);
            arb_mul(acb_imagref(e), ch, sn, p);
        }

        /* exp(alpha E - nu w - c) */
        acb_mul(t, ctx->alpha, e, p);
        acb_submul(t, ctx->nu, w, p);
        acb_sub(t, t, cc, p);
        acb_exp(t, t, p);

        arb_set_round(sh, wts + kk, p);
        acb_mul_arb(t, t, sh, p);
        acb_add(sum, sum, t, wp);
    }

    acb_mul(res, sum, h, wp);

    acb_clear(m); acb_clear(h); acb_clear(w); acb_clear(e);
    acb_clear(t); acb_clear(sum); acb_clear(cc);
    arb_clear(sh); arb_clear(ch); arb_clear(sn); arb_clear(cs);
}

/* the rule on the chord [a, b] between two vertices; lmax is an
   estimate of the maximum of Re(phi - c) on the chord */
static void
sdq_chord(acb_t res, sdq_ctx * ctx, cd c, cd a, cd b, int gi, double lmax)
{
    const void * x, * w;
    slong n = sdq_steps[gi];
    int N;

    if (ctx->N == 0)
    {
        acb_t aa, bb;
        acb_init(aa);
        acb_init(bb);
        acb_set_cd(aa, a);
        acb_set_cd(bb, b);
        sdq_chord_acb(res, ctx, c, aa, bb, a, b, gi);
        acb_clear(aa);
        acb_clear(bb);
        return;
    }

    N = (int) ceil((ctx->bits + fmin(lmax + 1.0, 0.0) / SD_LN2) / 53.0);
    N = FLINT_MAX(N, 1);
    N = FLINT_MIN(N, ctx->N);

    if (!ctx->par_set[N - 1])
    {
        if (N == 1) sdq_set_params_d1(&ctx->par[0].d1, ctx->alpha, ctx->nu);
        else if (N == 2) sdq_set_params_d2(&ctx->par[1].d2, ctx->alpha, ctx->nu);
        else if (N == 3) sdq_set_params_d3(&ctx->par[2].d3, ctx->alpha, ctx->nu);
        else sdq_set_params_d4(&ctx->par[3].d4, ctx->alpha, ctx->nu);
        ctx->par_set[N - 1] = 1;
    }

    sdq_get_nodes_dfloat(&x, &w, N, gi);
    if (N == 1) sdq_panel_d1(res, &ctx->par[0].d1, ctx->sym, c, a, b, x, w, n, ctx->scratch);
    else if (N == 2) sdq_panel_d2(res, &ctx->par[1].d2, ctx->sym, c, a, b, x, w, n, ctx->scratch);
    else if (N == 3) sdq_panel_d3(res, &ctx->par[2].d3, ctx->sym, c, a, b, x, w, n, ctx->scratch);
    else sdq_panel_d4(res, &ctx->par[3].d4, ctx->sym, c, a, b, x, w, n, ctx->scratch);
}

/* the chord from the starting point (a ball) to its approximation, the
   first vertex of the contour; the ellipses of all the chords from
   points of the ball are contained in disks around the vertex */
static int
sdq_start_chord(acb_t res, sdq_ctx * ctx, cd c, double logtol)
{
    const sd_problem * P = ctx->P;
    double ah = 0.5 * P->start_rad, logerr, nreq;
    acb_t b;
    mag_t err;
    int gi;

    acb_zero(res);
    if (ah == 0.0)
        return 1;

    gi = sdq_choose(&logerr, &nreq, ctx, P->start, ah, ah, c.re, logtol);
    if (gi < 0)
        return 0;

    acb_init(b);
    mag_init(err);
    acb_set_cd(b, P->start);
    sdq_chord_acb(res, ctx, c, ctx->start, b, P->start, P->start, gi);
    _sd_mag_set_exp_d(err, logerr);
    acb_add_error_mag(res, err);
    acb_clear(b);
    mag_clear(err);
    return 1;
}

/* add to res the integral over the segment [a, b], bisecting it until
   each part is negligible or a rule is found */
static int
sdq_bisect(acb_t res, sdq_ctx * ctx, cd c, cd a, cd b, double logtol, int depth)
{
    double g[2], mb, logerr, nreq;
    mag_t err;
    int gi, ok = 1;

    mag_init(err);
    mb = _sd_panel(g, ctx->P, a, b, c.re);

    if (mb <= logtol)
    {
        _sd_mag_set_exp_d(err, mb);
        acb_add_error_mag(res, err);
    }
    else if ((gi = sdq_choose_chord(&logerr, &nreq, ctx, a, b, c.re, logtol)) >= 0
        && (ctx->budget -= sdq_steps[gi]) >= 0)
    {
        cd phi, d1, d2;
        acb_t t;

        acb_init(t);
        _sd_phi(&phi, &d1, &d2, ctx->P, a);
        sdq_chord(t, ctx, c, a, b, gi, phi.re - c.re + 0.5 * cd_abs(d1) * cd_abs(cd_sub(b, a)));
        _sd_mag_set_exp_d(err, logerr);
        acb_add_error_mag(t, err);
        acb_add(res, res, t, ctx->wp);
        acb_clear(t);
    }
    else if (depth > 0 && ctx->budget >= 0)
    {
        cd m = cd_scale(cd_add(a, b), 0.5);
        ok = sdq_bisect(res, ctx, c, a, m, logtol - SD_LN2, depth - 1)
          && sdq_bisect(res, ctx, c, m, b, logtol - SD_LN2, depth - 1);
    }
    else
        ok = 0;

    mag_clear(err);
    return ok;
}

/* The quadrature along a piece: res (relative to exp(c)) with the tails
   and quadrature errors as radius; returns 0 if a chord cannot be
   resolved. */
static int
sdq_piece_integrate(acb_t res, sdq_ctx * ctx, const sd_piece * pc, double logtol, double cutoff)
{
    const sd_problem * P = ctx->P;
    const cd * pts = pc->pts;
    int npts = pc->npts, i, j, gi, gj, k, ok = 1;
    double * rephi, * mb, logerr, logerr2, nreq, nreq2, g[2], top = pc->c.re;
    mag_t err;
    acb_t t;

    acb_zero(res);
    mag_init(err);
    acb_init(t);

    for (k = 0; k < 2; k++)
    {
        _sd_mag_set_exp_d(err, pc->logtail[k]);
        acb_add_error_mag(res, err);
    }

    rephi = flint_malloc(sizeof(double) * 2 * (npts + 1));
    mb = rephi + npts + 1;

    /* the integrand at the vertices, and the majorants of the panels */
    for (i = 0; i < npts; i++)
    {
        cd phi, d1, d2;
        _sd_phi(&phi, &d1, &d2, P, pts[i]);
        rephi[i] = phi.re - top;
        if (i + 1 < npts)
            mb[i] = _sd_panel(g, P, pts[i], pts[i + 1], top);
    }

    /* greedy chords [i, j] */
    i = 0;
    while (i < npts - 1 && ok)
    {
        double lmax;

        /* panels below the tolerance are bounded, not integrated */
        if (mb[i] <= logtol)
        {
            double sum = mb[i];
            j = i + 1;
            while (j + 1 < npts && _sd_logaddexp(sum, mb[j]) <= logtol)
                sum = _sd_logaddexp(sum, mb[j++]);
            _sd_mag_set_exp_d(err, sum);
            acb_add_error_mag(res, err);
            i = j;
            continue;
        }

        j = i + 1;
        gi = sdq_choose_chord(&logerr, &nreq, ctx, pts[i], pts[j], top, logtol);
        if (gi < 0)
        {
            /* a long panel */
            ok = sdq_bisect(res, ctx, pc->c, pts[i], pts[j], logtol, 20);
            i = j;
            continue;
        }

        /* extend the chord over the following panels while this lowers
           its number of nodes per panel */
        while (j + 1 < npts && mb[j] > logtol)
        {
            /* the vertices between must stay close to the chord where
               the integrand is not negligible */
            for (k = i + 1; k <= j; k++)
                if (rephi[k] > -cutoff
                    && sdq_segdist(pts[k], pts[i], pts[j + 1]) > 0.3 * cd_abs(cd_sub(pts[k + 1], pts[k])))
                    break;
            if (k <= j)
                break;

            gj = sdq_choose_chord(&logerr2, &nreq2, ctx, pts[i], pts[j + 1], top, logtol);
            if (gj < 0 || nreq2 * (j - i) > nreq * (j + 1 - i))
                break;

            j++;
            gi = gj;
            logerr = logerr2;
            nreq = nreq2;
        }

        ctx->budget -= sdq_steps[gi];
        if (ctx->budget < 0)
        {
            ok = 0;
            break;
        }

        lmax = rephi[i];
        for (k = i + 1; k <= j; k++)
            lmax = fmax(lmax, rephi[k]);

        sdq_chord(t, ctx, pc->c, pts[i], pts[j], gi, lmax);
        _sd_mag_set_exp_d(err, logerr);
        acb_add_error_mag(t, err);
        acb_add(res, res, t, ctx->wp);

        i = j;
    }

    flint_free(rephi);
    mag_clear(err);
    acb_clear(t);
    return ok;
}

/* log of the radius of x */
static double
sdq_log_rad(const acb_t x)
{
    mag_t t;
    double r;

    mag_init(t);
    mag_hypot(t, arb_radref(acb_realref(x)), arb_radref(acb_imagref(x)));
    r = mag_get_d_log2_approx(t) * SD_LN2;
    mag_clear(t);
    return r;
}

int
_sd_integrate(acb_t res, const sd_contour * C, const sd_problem * P,
    const acb_t alpha, int sym, const acb_t nu, const acb_t start, slong prec)
{
    sdq_ctx ctx;
    cd cref;
    int i, k, ok = 1;
    double logL = C->logL[0], delta, cutoff, mag, wmax, xlo, xhi, qprec;
    acb_t sum, gsum, t, c;
    arf_t u;

    if (C->npieces == 0)
        return 0;

    /* The relative radius of the result due to the input balls will
       be about (r_alpha L1 + r_nu L2) / L0 in terms of the weighted
       L1 norms, plus the integrand at the starting point times its
       radius: more precision than this is wasted. */
    delta = _sd_logaddexp(sdq_log_rad(alpha) + C->logL[1], sdq_log_rad(nu) + C->logL[2]);
    if (start != NULL)
        delta = _sd_logaddexp(delta, sdq_log_rad(start) + C->piece[0].c.re);
    delta -= logL;
    qprec = fmin((double) prec, -delta / SD_LN2 + 10.0);
    qprec = fmax(qprec, 10.0);
    cutoff = fmin(qprec * SD_LN2 + 12.0, P->tailtol);

    /* the reference value, and the extent of the polygons where the
       integrand is not negligible */
    cref = C->piece[0].c;
    wmax = 0.0;
    xlo = INFINITY;
    xhi = -INFINITY;
    for (i = 0; i < C->npieces; i++)
    {
        const sd_piece * pc = C->piece + i;

        if (pc->c.re > cref.re)
            cref.re = pc->c.re;

        for (k = 0; k < pc->npts; k++)
        {
            cd phi, d1, d2;
            _sd_phi(&phi, &d1, &d2, P, pc->pts[k]);
            if (k == 0 || phi.re - pc->c.re > -cutoff)
            {
                wmax = fmax(wmax, cd_abs(pc->pts[k]));
                xlo = fmin(xlo, pc->pts[k].re);
                xhi = fmax(xhi, pc->pts[k].re);
            }
        }
    }

    /* the cancellation in phi(w) - c */
    mag = P->aal * exp(xhi) + ((P->abe != 0.0) ? P->abe * exp(-xlo) : 0.0)
        + P->anu * wmax + cd_abs(cref) + 2.0;
    if (!(mag < 1e300))
        return 0;

    ctx.P = P;
    ctx.alpha = alpha;
    ctx.nu = nu;
    ctx.start = start;
    ctx.sym = sym;
    ctx.bits = qprec + SDQ_GUARD + log2(mag);
    ctx.wp = (slong) ctx.bits + 8;
    ctx.par_set[0] = ctx.par_set[1] = ctx.par_set[2] = ctx.par_set[3] = 0;
    ctx.budget = 4000 + 40 * (slong) qprec;
    ctx.scratch = NULL;

    if (dfloat_is_supported() && ctx.bits <= 4 * 53)
    {
        ctx.N = FLINT_MAX(1, (int) ceil(ctx.bits / 53.0));
        ctx.wp = 53 * ctx.N + 10;
        ctx.imax = SDQ_DFLOAT_STEPS - 1;
        ctx.scratch = flint_malloc(SDQ_DFLOAT_NMAX * 12 * (ctx.N + 1) * sizeof(double));
    }
    else
    {
        ctx.N = 0;
        for (ctx.imax = SDQ_DFLOAT_STEPS - 1; ctx.imax + 1 < SDQ_STEPS
            && sdq_steps[ctx.imax + 1] <= 0.5 * qprec + 60.0; ctx.imax++) ;
    }

    acb_init(sum);
    acb_init(gsum);
    acb_init(t);
    acb_init(c);
    arf_init(u);

    /* sum the pieces: exp(c_piece - c_ref) times the piece integrals */
    for (i = 0; i < C->npieces && ok; i++)
    {
        const sd_piece * pc = C->piece + i;
        /* the tolerance for this piece, relative to exp(Re c) */
        double logtol = logL - (qprec + 12.0) * SD_LN2 - pc->c.re;

        ok = sdq_piece_integrate(t, &ctx, pc, logtol, cutoff);

        if (ok && i == 0 && start != NULL)
        {
            ok = sdq_start_chord(c, &ctx, pc->c, logtol);
            /* (the piece is traversed from the starting point: sign = 1) */
            acb_add(t, t, c, ctx.wp);
        }

        if (pc->sign < 0)
            acb_neg(t, t);

        arb_set_d(acb_realref(c), pc->c.re);
        arf_set_d(u, cref.re);
        arb_sub_arf(acb_realref(c), acb_realref(c), u, ctx.wp);
        arb_set_d(acb_imagref(c), pc->c.im);
        arf_set_d(u, cref.im);
        arb_sub_arf(acb_imagref(c), acb_imagref(c), u, ctx.wp);
        acb_exp(c, c, ctx.wp);
        acb_mul(t, t, c, ctx.wp);
        acb_add(pc->geom ? gsum : sum, pc->geom ? gsum : sum, t, ctx.wp);
    }

    /* the pieces with the factor 1 / (1 - q), q = exp(-2 pi i sigma nu) */
    if (ok && !acb_is_zero(gsum))
    {
        acb_const_pi(c, ctx.wp);
        acb_mul(c, c, nu, ctx.wp);
        acb_mul_onei(c, c);
        acb_mul_si(c, c, (P->lin_dir.im > 0.0) ? -2 : 2, ctx.wp);
        acb_exp(c, c, ctx.wp);
        acb_sub_ui(c, c, 1, ctx.wp);
        acb_neg(c, c);
        acb_div(gsum, gsum, c, ctx.wp);
        acb_add(sum, sum, gsum, ctx.wp);
    }

    if (ok)
    {
        acb_set_cd(c, cref);
        acb_exp(c, c, ctx.wp);
        acb_mul(res, sum, c, ctx.wp);

        ok = acb_is_finite(res);
    }

    flint_free(ctx.scratch);
    acb_clear(sum);
    acb_clear(gsum);
    acb_clear(t);
    acb_clear(c);
    arf_clear(u);
    return ok;
}
