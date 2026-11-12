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

    for all parameters in the input balls when phi and its derivatives
    are evaluated on the balls.

    What is rigorous and what is not.  The polygons, the choice of the
    chords, of the number of nodes and of X for each chord, the
    tolerances and the working precisions are heuristic: they are
    chosen in double precision from estimates (here and in sd_plan.c)
    so that the errors come out below the target, and the panels of the
    planned polygon are merged into longer chords while this does not
    increase the estimated number of nodes.  The result is then an
    enclosure whatever these choices are, because every term of its
    radius is computed rigorously for the choices made: the values of
    the rules in ball arithmetic with the input balls as parameters (so
    that the result encloses the integral on the whole input balls),
    the errors of the rules from the bound above with M computed by
    _sd_bound_disk, and the parts of the polygons where the integrand is
    estimated to be below the tolerance and the tails beyond the
    polygons, which are not integrated, by _sd_bound_segment and
    _sd_bound_tail (sd_bound.c, in arb and mag arithmetic).  The
    tolerance is relaxed to what the input balls permit.

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
    with dfloat, per node with acb), but not below what the cancellation
    at the node requires for the exponent to be accurate.
*/

#include "thread_support.h"
#include "arb_hypgeom.h"
#include "acb.h"
#include "acb_calc/impl.h"
#include "acb_hypgeom/sd.h"
#include "dfloat.h"

/* The rules used are those of the cache of acb_calc_integrate, with
   numbers of nodes in steps of about 2^(1/2) (steps of 2^(1/4) would
   save 8-10% of the nodes, at the price of about twice as many rules
   to precompute), up to about half the precision in bits, as in
   acb_calc_integrate. */
#define SDQ_STEPS ACB_CALC_GL_STEPS
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
    double cabs;             /* |c| of the reference value */
    double lref;             /* the level of Re(phi - c) below which the
                                precision is lowered: that of the largest
                                reference value, relative to the one of
                                the piece */
    slong wp;                /* working precision for acb arithmetic */
    int imax;                /* the largest rule */
    slong budget;            /* of nodes */
    sdq_params par[4];
    int par_set[4];
    void * scratch;
}
sdq_ctx;

/* Heuristic (in double precision; the error bound of the rule chosen
   is computed rigorously by sdq_gl_error).  The smallest rule (its
   index, or -1) whose error is estimated to be below exp(logtol) on any
   chord with half-length at most ah whose Bernstein ellipses with
   parameter X are contained in the disks |w - m| <= ah X + r0, with the
   best X, which is returned in *Xbest: the supremum of
   |exp(phi(w) - c)| (Re c = top) over the disk and the input balls is
   estimated as in the header comment.  Also sets *nreq to the number
   of nodes required if it could be any real number. */
static int
sdq_choose(double * Xbest, double * nreq, const sdq_ctx * ctx, cd m, double ah, double r0, double top, double logtol)
{
    static const double Xs[] = { 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0,
        1e2, 1e3, 1e4, 1e6, 1e8, 1e10, 1e12 };
    const sd_problem * P = ctx->P;
    cd phi, d1, d2;
    double am = cd_abs(m), A, B, D, base, bestn = INFINITY;
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
        *Xbest = X;
    }

    *nreq = fmax(bestn, 1.0);

    for (i = 0; i <= ctx->imax; i++)
        if (sdq_steps[i] >= bestn)
            break;
    if (i > ctx->imax)
        return -1;

    return i;
}

static int
sdq_choose_chord(double * Xbest, double * nreq, const sdq_ctx * ctx, cd a, cd b, double top, double logtol)
{
    return sdq_choose(Xbest, nreq, ctx, cd_scale(cd_add(a, b), 0.5),
        0.5 * cd_abs(cd_sub(b, a)), 0.0, top, logtol);
}

/* Rigorous.  err = a bound for the error of the n-point Gauss-Legendre
   rule for the integral of exp(phi(w) - c) over any chord with
   half-length at most ah whose midpoint is within r0 of m (a ball),
   for all parameters in the balls:

       (64 / 15) M ah / ((rho - 1) rho^(2 n - 1)),   rho = X + sqrt(X^2 - 1),

   where M bounds the integrand on the disk |w - m| <= ah X + r0, which
   contains the Bernstein ellipses with parameter rho of these chords
   (X > 1 is an exact number). */
static void
sdq_gl_error(mag_t err, const sdq_ctx * ctx, const acb_t c, const acb_t m,
    const mag_t ah, const mag_t r0, double X, slong n)
{
    arb_t rho, t;
    mag_t R, M;

    arb_init(rho);
    arb_init(t);
    mag_init(R);
    mag_init(M);

    if (!(X > 1.0) || !(X < 1e300))
    {
        mag_inf(err);
    }
    else
    {
        /* R = ah X + r0 */
        mag_set_d(R, X);
        mag_mul(R, R, ah);
        mag_add(R, R, r0);
        _sd_bound_disk(M, ctx->alpha, ctx->sym, ctx->nu, m, R, c);

        /* (rho - 1) rho^(2 n - 1) */
        arb_set_d(rho, X);
        arb_mul(t, rho, rho, 64);
        arb_sub_ui(t, t, 1, 64);
        arb_sqrt(t, t, 64);
        arb_add(rho, rho, t, 64);
        arb_pow_ui(t, rho, 2 * n - 1, 64);
        arb_sub_ui(rho, rho, 1, 64);
        arb_mul(t, t, rho, 64);
        arb_ui_div(t, 64, t, 64);
        arb_div_ui(t, t, 15, 64);

        arb_get_mag(err, t);
        mag_mul(err, err, M);
        mag_mul(err, err, ah);
    }

    arb_clear(rho);
    arb_clear(t);
    mag_clear(R);
    mag_clear(M);
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
   evaluated with the corresponding number of bits less.  The nodes are
   evaluated in parallel when there are several threads and enough
   work for them. */
typedef struct
{
    const sdq_ctx * ctx;
    cd c, md, hd;
    acb_srcptr m, h, cc;
    arb_srcptr x, wts;
    slong n, wp;
    acb_ptr v;
}
sdq_acb_work;

/* The number of bits lost to cancellation in phi(w) - c at w, plus 10:
   the smallest precision at which the exponent is accurate enough for
   the exponential not to inflate its error.  (Where the integrand is
   small, the precision is lowered by its size; far out in the valleys
   the terms of phi are larger than where the integrand is not
   negligible, which is what the working precision allows for.) */
static double
sdq_min_bits(const sdq_ctx * ctx, cd w)
{
    const sd_problem * P = ctx->P;
    double sz = P->aal * exp(w.re) + ((P->abe != 0.0) ? P->abe * exp(-w.re) : 0.0)
        + P->anu * cd_abs(w) + ctx->cabs + 2.0;

    if (!(sz < 1e300))
        return 1100.0;
    return log2(sz) + 10.0;
}

/* t = the term of the node k; w, e are temporaries */
static void
sdq_node_acb(acb_t t, const sdq_acb_work * W, slong k, acb_t w, acb_t e, arb_ptr r)
{
    const sdq_ctx * ctx = W->ctx;
    slong n = W->n, kk, p;
    cd phi, d1, d2, wd;
    double xd;
    arb_ptr sh = r, ch = r + 1, sn = r + 2, cs = r + 3;

    kk = (2 * k < n) ? k : n - 1 - k;
    xd = arf_get_d(arb_midref(W->x + kk), ARF_RND_NEAR);
    if (2 * k >= n)
        xd = -xd;

    wd = cd_add(W->md, cd_scale(W->hd, xd));
    _sd_phi(&phi, &d1, &d2, ctx->P, wd);
    p = W->wp + (slong) fmin((phi.re - W->c.re - ctx->lref) / SD_LN2 + 2.0, 0.0);
    p = FLINT_MAX(p, 32);
    p = FLINT_MAX(p, (slong) sdq_min_bits(ctx, wd));

    /* w = m +- h x (the cached nodes can have a higher precision) */
    arb_set_round(sh, W->x + kk, p);
    acb_mul_arb(w, W->h, sh, p);
    if (2 * k >= n)
        acb_sub(w, W->m, w, p);
    else
        acb_add(w, W->m, w, p);

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
    acb_sub(t, t, W->cc, p);
    acb_exp(t, t, p);

    arb_set_round(sh, W->wts + kk, p);
    acb_mul_arb(t, t, sh, p);
}

static void
sdq_node_worker(slong k, void * arg)
{
    const sdq_acb_work * W = (const sdq_acb_work *) arg;
    acb_t w, e;
    arb_struct r[4];
    int i;

    acb_init(w); acb_init(e);
    for (i = 0; i < 4; i++)
        arb_init(r + i);

    sdq_node_acb(W->v + k, W, k, w, e, r);

    acb_clear(w); acb_clear(e);
    for (i = 0; i < 4; i++)
        arb_clear(r + i);
}

static void
sdq_chord_acb(acb_t res, const sdq_ctx * ctx, cd c, const acb_t a, const acb_t b,
    cd ad, cd bd, int gi)
{
    sdq_acb_work W;
    slong n = sdq_steps[gi], k, wp = ctx->wp;
    acb_t m, h, w, e, t, sum, cc;
    arb_struct r[4];
    int i;

    acb_init(m); acb_init(h); acb_init(w); acb_init(e);
    acb_init(t); acb_init(sum); acb_init(cc);
    for (i = 0; i < 4; i++)
        arb_init(r + i);

    _acb_calc_gl_nodes(&W.x, &W.wts, gi, wp);

    acb_add(m, a, b, wp);
    acb_mul_2exp_si(m, m, -1);
    acb_sub(h, b, a, wp);
    acb_mul_2exp_si(h, h, -1);
    acb_set_cd(cc, c);

    W.ctx = ctx;
    W.c = c;
    W.md = cd_scale(cd_add(ad, bd), 0.5);
    W.hd = cd_scale(cd_sub(bd, ad), 0.5);
    W.m = m; W.h = h; W.cc = cc;
    W.n = n;
    W.wp = wp;

    /* in parallel if the work is well above the overhead of starting
       the threads (some microseconds per node at a few hundred bits) */
    if (flint_get_num_threads() >= 2 && n >= 8 && (double) n * (double) wp >= 30000.0)
    {
        W.v = _acb_vec_init(n);
        flint_parallel_do(sdq_node_worker, &W, n, -1, FLINT_PARALLEL_STRIDED);
        /* in the order of the serial loop: the result does not depend
           on the number of threads */
        for (k = 0; k < n; k++)
            acb_add(sum, sum, W.v + k, wp);
        _acb_vec_clear(W.v, n);
    }
    else
    {
        for (k = 0; k < n; k++)
        {
            sdq_node_acb(t, &W, k, w, e, r);
            acb_add(sum, sum, t, wp);
        }
    }

    acb_mul(res, sum, h, wp);

    acb_clear(m); acb_clear(h); acb_clear(w); acb_clear(e);
    acb_clear(t); acb_clear(sum); acb_clear(cc);
    for (i = 0; i < 4; i++)
        arb_clear(r + i);
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

    N = (int) ceil(fmax(ctx->bits + fmin(lmax - ctx->lref + 1.0, 0.0) / SD_LN2,
            fmax(sdq_min_bits(ctx, a), sdq_min_bits(ctx, b))) / 53.0);
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
   first vertex of the contour; the midpoints of all the chords from
   points of the ball are within their half-length of the vertex */
static int
sdq_start_chord(acb_t res, sdq_ctx * ctx, cd c, double logtol)
{
    const sd_problem * P = ctx->P;
    double X, nreq;
    acb_t b, cc;
    mag_t err, ah;
    int gi, ok = 0;

    acb_zero(res);

    acb_init(b);
    acb_init(cc);
    mag_init(err);
    mag_init(ah);

    /* the largest half-length */
    acb_set_cd(b, P->start);
    acb_sub(cc, ctx->start, b, 64);
    acb_get_mag(ah, cc);
    mag_mul_2exp_si(ah, ah, -1);

    if (mag_is_zero(ah))
    {
        ok = 1;
    }
    else if (mag_is_finite(ah))
    {
        double ahd = mag_get_d(ah);

        gi = sdq_choose(&X, &nreq, ctx, P->start, ahd, ahd, c.re, logtol);
        if (gi >= 0)
        {
            sdq_chord_acb(res, ctx, c, ctx->start, b, P->start, P->start, gi);
            acb_set_cd(cc, c);
            sdq_gl_error(err, ctx, cc, b, ah, ah, X, sdq_steps[gi]);
            acb_add_error_mag(res, err);
            ok = 1;
        }
    }

    acb_clear(b);
    acb_clear(cc);
    mag_clear(err);
    mag_clear(ah);
    return ok;
}

/* The quadrature along a piece: res (relative to exp(c)) with the tails
   and quadrature errors as radius; returns 0 if a tail is not certified
   or no rule is found for a chord.  The vertices where the tails begin
   and the chords are chosen heuristically; the bounds for the tails and
   for the errors of the rules are rigorous. */
static int
sdq_piece_integrate(acb_t res, sdq_ctx * ctx, const sd_piece * pc, double logtol, double cutoff)
{
    const sd_problem * P = ctx->P;
    const cd * pts = pc->pts;
    int npts = pc->npts, i, j, gi, gj, k, ok = 1, i0, i1;
    double * rephi, * mb, X, X2, nreq, nreq2, g[2], top = pc->c.re, sum;
    mag_t err, ah, zero;
    acb_t t, c, m;

    acb_zero(res);
    mag_init(err);
    mag_init(ah);
    mag_init(zero);
    acb_init(t);
    acb_init(c);
    acb_init(m);

    acb_set_cd(c, pc->c);

    rephi = flint_malloc(sizeof(double) * 2 * (npts + 1));
    mb = rephi + npts + 1;

    /* estimates of the integrand at the vertices and of the integrals
       over the panels */
    for (i = 0; i < npts; i++)
    {
        cd phi, d1, d2;
        _sd_phi(&phi, &d1, &d2, P, pts[i]);
        rephi[i] = phi.re - top;
        if (i + 1 < npts)
            mb[i] = _sd_panel(g, P, pts[i], pts[i + 1], top);
    }

    /* The ends: the polygon is integrated between the outermost
       vertices i0 and i1 beyond which the integral is estimated to be
       below the tolerance; the rest, with the tails into the valleys,
       is bounded. */
    i0 = 0;
    i1 = npts - 1;
    /* (logtail is -inf at an end which is not in a valley) */
    sum = pc->logtail[0];
    while (i0 + 1 < npts && _sd_logaddexp(sum, mb[i0]) <= logtol)
        sum = _sd_logaddexp(sum, mb[i0++]);
    sum = pc->logtail[1];
    while (i1 > i0 && _sd_logaddexp(sum, mb[i1 - 1]) <= logtol)
        sum = _sd_logaddexp(sum, mb[--i1]);
    for (k = 0; k < 2 && ok; k++)
    {
        ok = _sd_bound_end(err, ctx->alpha, ctx->sym, ctx->nu, pc, c, k ? i1 : i0, k, P->lin_dir, P, logtol - 2.0);
        acb_add_error_mag(res, err);
    }

    /* greedy chords [i, j] */
    i = i0;
    while (i < i1 && ok)
    {
        double lmax;

        j = i + 1;

        /* a panel estimated to be below the tolerance is bounded, not
           integrated */
        if (mb[i] <= logtol)
        {
            acb_set_cd(m, pts[i]);
            acb_set_cd(t, pts[j]);
            _sd_bound_segment(err, 1, ctx->alpha, ctx->sym, ctx->nu, m, t, c);
            acb_add_error_mag(res, err);
            i = j;
            continue;
        }

        gi = sdq_choose_chord(&X, &nreq, ctx, pts[i], pts[j], top, logtol);
        if (gi < 0)
        {
            /* No rule is found for this panel.  This only happens
               for input balls so wide that the result would have no
               accuracy; the integration fails. */
            ok = 0;
            break;
        }

        /* extend the chord over the following panels while this lowers
           its number of nodes per panel */
        while (j < i1 && mb[j] > logtol)
        {
            /* the vertices between must stay close to the chord where
               the integrand is not negligible */
            for (k = i + 1; k <= j; k++)
                if (rephi[k] > -cutoff
                    && sdq_segdist(pts[k], pts[i], pts[j + 1]) > 0.3 * cd_abs(cd_sub(pts[k + 1], pts[k])))
                    break;
            if (k <= j)
                break;

            gj = sdq_choose_chord(&X2, &nreq2, ctx, pts[i], pts[j + 1], top, logtol);
            if (gj < 0 || nreq2 * (j - i) > nreq * (j + 1 - i))
                break;

            j++;
            gi = gj;
            X = X2;
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

        /* the error of the rule: the midpoint and the half-length of
           the chord between two doubles */
        acb_set_cd(m, pts[i]);
        {
            acb_t b;
            acb_init(b);
            acb_set_cd(b, pts[j]);
            acb_sub(b, b, m, 128);
            acb_mul_2exp_si(b, b, -1);
            acb_get_mag(ah, b);
            acb_add(m, m, b, 128);
            acb_clear(b);
        }
        sdq_gl_error(err, ctx, c, m, ah, zero, X, sdq_steps[gi]);
        acb_add_error_mag(t, err);
        acb_add(res, res, t, ctx->wp);

        i = j;
    }

    flint_free(rephi);
    mag_clear(err);
    mag_clear(ah);
    mag_clear(zero);
    acb_clear(t);
    acb_clear(c);
    acb_clear(m);
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
    return _sd_integrate_pieces(res, C, P, alpha, sym, nu, start, prec, NULL);
}

int
_sd_integrate_pieces(acb_t res, const sd_contour * C, const sd_problem * P,
    const acb_t alpha, int sym, const acb_t nu, const acb_t start, slong prec, const char * sel)
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
    /* The errors are to be small compared to the integral, which is
       about exp(logL): narrower than exp(Re c_ref) by the width of the
       saddle. */
    ctx.bits = qprec + SDQ_GUARD + log2(mag) + fmax((cref.re - logL) / SD_LN2, 0.0);
    ctx.cabs = cd_abs(cref);
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

        if (sel != NULL && !sel[i])
            continue;

        ctx.lref = cref.re - pc->c.re;

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
        /* (the geometric series must converge: |q| < 1) */
        {
            mag_t qm;
            mag_init(qm);
            acb_get_mag(qm, c);
            ok = mag_cmp_2exp_si(qm, 0) < 0;
            mag_clear(qm);
        }
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
