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

    For large parameters, the vertices are doubles in the frame of a
    saddle (see sd.h), and the integrand is evaluated in the frame, by
    a formula without that cancellation (sdq_set_frame); the error
    bounds are computed from the exact points and reference values as
    before, at a higher precision.
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
#define SDQ_DFLOAT_STEPS ACB_CALC_GL_DFLOAT_STEPS
#define SDQ_DFLOAT_NMAX ACB_CALC_GL_DFLOAT_NMAX
/* bits of guard beyond prec and the cancellation */
#define SDQ_GUARD 4.0
/* more bits for the points in a frame with rel = 1 */
#define SDQ_REL_GUARD 16
/* the largest number of Taylor coefficients of the exponent with dfloat */
#define SDQ_REL_NC 24

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
    const sd_frame * F;      /* the frame of the piece being integrated */
    /* for a frame with rel = 1 (see sdq_set_frame): */
    acb_struct ka, kb, kr, kk;  /* a, b, r and phi(w0) - cx for the input balls */
    acb_struct ks, kd;       /* a + b, a - b */
    acb_ptr coef;            /* the Taylor coefficients of the exponent in
                                the frame (from the second), for dfloat */
    acb_struct shift;        /* the reference value of the piece minus cx */
    const sd_frame * kF;     /* the frame these are for */
    slong hp;                /* the precision for computing them */
    int N;                   /* the largest dfloat precision, or 0 for acb */
    double bits;             /* prec and guard bits, plus the cancellation */
    double cabs;             /* |c| of the reference value (in the frame) */
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
    const sd_frame * F = ctx->F;
    cd phi, d1, d2;
    double am = cd_abs(cd_add(F->w0, m)), A, B, D, base, bestn = INFINITY;
    int ix, i;

    _sd_phi(&phi, &d1, &d2, P, F, m);
    A = cd_abs(d2);
    B = cd_abs(_sd_phi3(P, F, m, d1));
    D = cd_abs(d1);
    base = phi.re - top + log(64.0 / 15.0) + log(ah) - logtol + P->rnu * am;
    /* (the rounding of the estimate) */
    if (!F->rel)
        base += 64.0 * SD_EPS * (P->anu * am + fabs(top) + 1.0);
    else
        base += 4.0 * _sd_margin(P, F, m) + 64.0 * SD_EPS * (fabs(top) + 1.0);

    for (ix = 0; ix < (int) (sizeof(Xs) / sizeof(double)); ix++)
    {
        double X = Xs[ix], rho = X + sqrt(X * X - 1.0), lrho = log(rho), R = ah * X + r0;
        double sh, chm1, shmr, need, n, size, ball, round;

        if (R > 300.0)
            break;

        /* cosh R - 1 and sinh R - R, rounded up */
        sh = sinh(0.5 * R);
        chm1 = 2.0 * sh * sh;
        shmr = (R < 0.5) ? R * R * R * 0.17 : sinh(R) - R;
        size = _sd_size(F, m.re - R, m.re + R);
        ball = _sd_ball(F, m.re - R, m.re + R);
        /* (the rounding of the estimate: of the terms of phi, or in a
           frame with rel = 1 of those of its variation) */
        if (!F->rel)
            round = 64.0 * SD_EPS * (size + P->anu * R);
        else
            round = 64.0 * SD_EPS * (D * R + A * chm1 + B * shmr);

        /* (2n-1) log rho >= log(64/15) + log M + log(ah) - log(rho-1) - logtol */
        need = base + (D * R + A * chm1 + B * shmr) * (1.0 + 1e-12)
            + ball + P->rnu * R + round - log(rho - 1.0);
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
    mag_t R, M;

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
        _sd_bound_disk(M, ctx->alpha, ctx->sym, ctx->nu, m, R, c, ctx->P->bprec);
        _acb_calc_gl_error_bound(err, M, ah, X, n);
    }

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

/* The frame of the pieces to be integrated next, and their reference
   value.  For a frame with rel = 1, the integrand is evaluated at the
   point w0 + w of the frame as

       exp(a (e^w - 1 - w) + b (e^-w - 1 + w) + r w + (phi(w0) - cx) - shift)

   (times exp(cx + shift), where cx + shift is the reference value of
   the piece), with balls for a = alpha e^w0, b = beta e^-w0,
   r = a - b - nu and phi(w0) - cx = a + b - nu w0 - cx computed here from
   the input balls, at a precision at which the cancellation in the last
   two leaves the working precision.  The terms of the exponent are
   then of the size of the exponent itself near the saddle, whatever
   the size of the parameters, and so is the working precision.  For
   small w, where e^w - 1 - w cancels (and its evaluation on a ball
   for w would lose even more), the first two terms are evaluated as

       (a + b) (cosh w - 1) + (a - b) (sinh w - w)

   by the Taylor series (see sdq_node_acb). */
static void
sdq_set_frame(sdq_ctx * ctx, const sd_piece * pc)
{
    const sd_frame * F = pc->fr;
    slong hp = ctx->hp;
    acb_t e;

    ctx->F = F;
    if (!F->rel)
        return;

    acb_sub(&ctx->shift, pc->cx, &F->cx, ARF_PREC_EXACT);

    if (ctx->kF == F)
        return;

    acb_init(e);
    acb_exp(e, &F->base, hp);
    acb_mul(&ctx->ka, ctx->alpha, e, hp);
    acb_zero(&ctx->kb);
    if (ctx->sym)
    {
        acb_div(&ctx->kb, ctx->alpha, e, hp);
        acb_neg(&ctx->kb, &ctx->kb);
    }
    acb_add(&ctx->ks, &ctx->ka, &ctx->kb, hp);
    acb_sub(&ctx->kd, &ctx->ka, &ctx->kb, hp);
    acb_sub(&ctx->kr, &ctx->kd, ctx->nu, hp);
    acb_set(&ctx->kk, &ctx->ks);
    acb_submul(&ctx->kk, ctx->nu, &F->base, hp);
    acb_sub(&ctx->kk, &ctx->kk, &F->cx, hp);
    ctx->kF = F;
    acb_clear(e);

    /* r w + (a + b) (cosh w - 1) + (a - b) (sinh w - w) */
    if (ctx->N != 0)
    {
        slong k;
        arb_t f;

        arb_init(f);
        arb_one(f);
        acb_set(ctx->coef + 1, &ctx->kr);
        for (k = 2; k < SDQ_REL_NC; k++)
        {
            arb_div_ui(f, f, k, hp);
            acb_mul_arb(ctx->coef + k, (k % 2) ? &ctx->kd : &ctx->ks, f, hp);
        }
        arb_clear(f);
    }
}

/* the vertex j of a piece in its frame (exact) */
static void
sdq_vertex(acb_t w, const sd_piece * pc, slong j)
{
    if (j == 0 && pc->v0 != NULL)
        acb_sub(w, pc->v0, &pc->fr->base, ARF_PREC_EXACT);
    else
        acb_set_d_d(w, pc->pts[j].re, pc->pts[j].im);
}

/* The rule sdq_steps[gi] on the chord [a, b] in acb arithmetic.  The
   endpoints are balls (points of the frame); ad and bd are their double
   approximations, used
   to estimate the size of the integrand at each node, which is then
   evaluated with the corresponding number of bits less.  The nodes are
   evaluated in parallel when there are several threads and enough
   work for them. */
typedef struct
{
    const sdq_ctx * ctx;
    double top;
    cd md, hd;
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
    const sd_frame * F = ctx->F;
    double sz = P->aal * exp(w.re) + ((P->abe != 0.0) ? P->abe * exp(-w.re) : 0.0)
        + P->anu * cd_abs(w) + ctx->cabs + 2.0;

    /* (the terms of the exponent in the frame) */
    if (F->rel)
        sz = _sd_relsize(F, w) + ctx->cabs + 2.0;

    if (!(sz < 1e300))
        return 1100.0;
    return log2(sz) + 10.0;
}

/* t = the term of the node k; w, e are temporaries */
static void
sdq_node_acb(acb_t t, const sdq_acb_work * W, slong k, acb_t w, acb_t e, arb_ptr r)
{
    const sdq_ctx * ctx = W->ctx;
    slong n = W->n, kk, p, q;
    cd phi, d1, d2, wd;
    double xd;
    arb_ptr sh = r, ch = r + 1, sn = r + 2, cs = r + 3;

    kk = (2 * k < n) ? k : n - 1 - k;
    xd = arf_get_d(arb_midref(W->x + kk), ARF_RND_NEAR);
    if (2 * k >= n)
        xd = -xd;

    wd = cd_add(W->md, cd_scale(W->hd, xd));
    _sd_phi(&phi, &d1, &d2, ctx->P, ctx->F, wd);
    p = W->wp + (slong) fmax(fmin((phi.re - W->top - ctx->lref) / SD_LN2 + 2.0, 0.0), -1e9);
    p = FLINT_MAX(p, 32);
    q = (slong) sdq_min_bits(ctx, wd);
    p = FLINT_MAX(p, q);

    /* w = m +- h x (the cached nodes can have a higher precision) */
    if (ctx->F->rel)
    {
        const sd_frame * F = ctx->F;
        double mw = cd_abs(wd);
        slong pw = p + SDQ_REL_GUARD;
        acb_t u;

        /* (m and h are exact for a chord between two vertices) */
        arb_set_round(sh, W->x + kk, pw);
        acb_mul_arb(w, W->h, sh, pw);
        if (2 * k >= n)
            acb_sub(w, W->m, w, pw);
        else
            acb_add(w, W->m, w, pw);

        acb_init(u);

        if (mw < 0.125)
        {
            /* With z = w^2: cosh w - 1 = sum_{j >= 1} z^j / (2 j)! and
               sinh w - w = w sum_{j >= 1} z^j / (2 j + 1)!, with J
               terms; the rest of each sum is at most twice its first
               term, at most 2 |z|^(J + 1) / (2 J + 2)!.  J is chosen
               so that this is negligible for the exponent. */
            double z2 = mw * mw, tb = (F->aa + F->ab) * z2 * z2 / 12.0;
            slong J = 1, j;
            mag_t zb, fb;

            while (tb > ldexp(1.0, -(int) FLINT_MIN(p + 4, 1000)) && J < 500)
            {
                J++;
                tb *= z2 / ((2.0 * J + 1.0) * (2.0 * J + 2.0));
            }

            mag_init(zb);
            mag_init(fb);

            acb_mul(e, w, w, pw);
            acb_zero(t);
            acb_zero(u);
            for (j = J; j >= 1; j--)
            {
                acb_add_ui(t, t, 1, pw);
                acb_mul(t, t, e, pw);
                acb_div_ui(t, t, (2 * j - 1) * (2 * j), pw);
                acb_add_ui(u, u, 1, pw);
                acb_mul(u, u, e, pw);
                acb_div_ui(u, u, (2 * j) * (2 * j + 1), pw);
            }

            acb_get_mag(zb, e);
            mag_pow_ui(zb, zb, J + 1);
            mag_rfac_ui(fb, 2 * J + 2);
            mag_mul(zb, zb, fb);
            mag_mul_2exp_si(zb, zb, 1);
            acb_add_error_mag(t, zb);
            acb_add_error_mag(u, zb);

            acb_mul(u, u, w, pw);
            acb_mul(t, t, &ctx->ks, p);
            acb_addmul(t, u, &ctx->kd, p);

            mag_clear(zb);
            mag_clear(fb);
        }
        else
        {
            /* a (e^w - 1 - w) + b (e^-w - 1 + w), where
               e^-w - 1 = -(e^w - 1) / e^w; the subtractions lose at
               most a few bits */
            acb_expm1(e, w, pw);
            acb_sub(t, e, w, pw);
            acb_mul(t, t, &ctx->ka, p);
            if (ctx->sym)
            {
                acb_add_ui(u, e, 1, pw);
                acb_div(u, e, u, pw);
                acb_sub(u, w, u, pw);
                acb_addmul(t, u, &ctx->kb, p);
            }
        }

        acb_addmul(t, &ctx->kr, w, p);
        acb_add(t, t, &ctx->kk, p);
        acb_sub(t, t, W->cc, p);
        acb_exp(t, t, p);

        arb_set_round(sh, W->wts + kk, p);
        acb_mul_arb(t, t, sh, p);
        acb_clear(u);
        return;
    }

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
sdq_chord_acb(acb_t res, const sdq_ctx * ctx, const acb_t cc, double top, const acb_t a, const acb_t b,
    cd ad, cd bd, int gi)
{
    sdq_acb_work W;
    slong n = sdq_steps[gi], k, wp = ctx->wp;
    acb_t m, h, w, e, t, sum;
    arb_struct r[4];
    int i;

    acb_init(m); acb_init(h); acb_init(w); acb_init(e);
    acb_init(t); acb_init(sum);
    for (i = 0; i < 4; i++)
        arb_init(r + i);

    /* (in a frame with rel = 1, the points are needed with a few more
       bits; the vertices are small numbers, and their sums are exact) */
    if (ctx->F->rel)
    {
        slong pm = (acb_is_exact(a) && acb_is_exact(b)) ? ARF_PREC_EXACT : wp + SDQ_REL_GUARD;

        _acb_calc_gl_nodes(&W.x, &W.wts, gi, wp + SDQ_REL_GUARD);
        acb_add(m, a, b, pm);
        acb_mul_2exp_si(m, m, -1);
        acb_sub(h, b, a, pm);
        acb_mul_2exp_si(h, h, -1);
    }
    else
    {
        _acb_calc_gl_nodes(&W.x, &W.wts, gi, wp);

        acb_add(m, a, b, wp);
        acb_mul_2exp_si(m, m, -1);
        acb_sub(h, b, a, wp);
        acb_mul_2exp_si(h, h, -1);
    }

    W.ctx = ctx;
    W.top = top;
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
    acb_clear(t); acb_clear(sum);
    for (i = 0; i < 4; i++)
        arb_clear(r + i);
}

/* the value to subtract from the exponent in sdq_chord_acb, for the
   pieces with the reference value c */
static void
sdq_get_cc(acb_t cc, const sdq_ctx * ctx, cd c)
{
    if (ctx->F->rel)
        acb_set(cc, &ctx->shift);
    else
        acb_set_cd(cc, c);
}

/* the rule on the chord between the vertices i and j of the piece;
   lmax is an estimate of the maximum of Re(phi - c) on the chord */
static void
sdq_chord(acb_t res, sdq_ctx * ctx, const sd_piece * pc, int i, int j, int gi, double lmax)
{
    const void * x, * w;
    slong n = sdq_steps[gi];
    cd c = pc->c, a = pc->pts[i], b = pc->pts[j];
    int N;

    if (ctx->N == 0)
    {
        acb_t aa, bb, cc;
        acb_init(aa);
        acb_init(bb);
        acb_init(cc);
        sdq_vertex(aa, pc, i);
        sdq_vertex(bb, pc, j);
        sdq_get_cc(cc, ctx, c);
        sdq_chord_acb(res, ctx, cc, pc->top, aa, bb, a, b, gi);
        acb_clear(aa);
        acb_clear(bb);
        acb_clear(cc);
        return;
    }

    N = (int) ceil(fmax(ctx->bits + fmin(lmax - ctx->lref + 1.0, 0.0) / SD_LN2,
            fmax(sdq_min_bits(ctx, a), sdq_min_bits(ctx, b))) / 53.0);
    N = FLINT_MAX(N, 1);
    N = FLINT_MIN(N, ctx->N);

    if (ctx->F->rel)
    {
        /* The Taylor polynomial of the exponent, with J terms of each
           of the series of cosh w - 1 and sinh w - w in w^2 (see
           sdq_node_acb): for all points of the chord, |w| <= wmax, and
           the rest is at most
           (|a + b| + |a - b| |w|) 2 |w|^(2 J + 2) / (2 J + 2)!. */
        const sd_frame * F = ctx->F;
        double wmax = fmax(cd_abs(a), cd_abs(b)), z2 = wmax * wmax;
        double tb = (F->aa + F->ab) * z2 * z2 / 12.0;
        slong J = 1;

        while (tb > ldexp(1.0, -(53 * N + 4)) && J < SDQ_REL_NC)
        {
            J++;
            tb *= z2 / ((2.0 * J + 1.0) * (2.0 * J + 2.0));
        }

        if (wmax < 0.125 && 2 * J + 2 <= SDQ_REL_NC && !(i == 0 && pc->v0 != NULL))
        {
            mag_t t, u, v;
            acb_t c0;

            mag_init(t);
            mag_init(u);
            mag_init(v);
            acb_init(c0);

            mag_set_d(t, wmax);
            acb_get_mag(u, &ctx->kd);
            mag_mul(u, u, t);
            acb_get_mag(v, &ctx->ks);
            mag_add(u, u, v);
            mag_pow_ui(t, t, 2 * J + 2);
            mag_mul(u, u, t);
            mag_rfac_ui(t, 2 * J + 2);
            mag_mul(u, u, t);
            mag_mul_2exp_si(u, u, 1);

            acb_sub(c0, &ctx->kk, &ctx->shift, ctx->hp);
            acb_add_error_mag(c0, u);
            acb_swap(ctx->coef, c0);

            _acb_calc_gl_exp_poly_dfloat(res, ctx->coef, 2 * J + 2, a.re, a.im, b.re, b.im,
                gi, N, ctx->scratch);

            mag_clear(t);
            mag_clear(u);
            mag_clear(v);
            acb_clear(c0);
        }
        else
        {
            acb_t aa, bb, cc;
            acb_init(aa);
            acb_init(bb);
            acb_init(cc);
            sdq_vertex(aa, pc, i);
            sdq_vertex(bb, pc, j);
            sdq_get_cc(cc, ctx, c);
            sdq_chord_acb(res, ctx, cc, pc->top, aa, bb, a, b, gi);
            acb_clear(aa);
            acb_clear(bb);
            acb_clear(cc);
        }
        return;
    }

    if (!ctx->par_set[N - 1])
    {
        if (N == 1) sdq_set_params_d1(&ctx->par[0].d1, ctx->alpha, ctx->nu);
        else if (N == 2) sdq_set_params_d2(&ctx->par[1].d2, ctx->alpha, ctx->nu);
        else if (N == 3) sdq_set_params_d3(&ctx->par[2].d3, ctx->alpha, ctx->nu);
        else sdq_set_params_d4(&ctx->par[3].d4, ctx->alpha, ctx->nu);
        ctx->par_set[N - 1] = 1;
    }

    _acb_calc_gl_nodes_dfloat(&x, &w, N, gi);
    if (N == 1) sdq_panel_d1(res, &ctx->par[0].d1, ctx->sym, c, a, b, x, w, n, ctx->scratch);
    else if (N == 2) sdq_panel_d2(res, &ctx->par[1].d2, ctx->sym, c, a, b, x, w, n, ctx->scratch);
    else if (N == 3) sdq_panel_d3(res, &ctx->par[2].d3, ctx->sym, c, a, b, x, w, n, ctx->scratch);
    else sdq_panel_d4(res, &ctx->par[3].d4, ctx->sym, c, a, b, x, w, n, ctx->scratch);
}

/* the chord from the starting point (a ball) to its approximation, the
   first vertex of the contour; the midpoints of all the chords from
   points of the ball are within their half-length of the vertex */
static int
sdq_start_chord(acb_t res, sdq_ctx * ctx, const sd_piece * pc, double logtol)
{
    double X, nreq;
    acb_t a, b, cc;
    mag_t err, ah;
    int gi, ok = 0;

    acb_zero(res);

    acb_init(a);
    acb_init(b);
    acb_init(cc);
    mag_init(err);
    mag_init(ah);

    /* the largest half-length */
    _sd_piece_vertex(b, pc, 0);
    acb_sub(cc, ctx->start, b, ctx->P->bprec);
    acb_get_mag(ah, cc);
    mag_mul_2exp_si(ah, ah, -1);

    if (mag_is_zero(ah))
    {
        ok = 1;
    }
    else if (mag_is_finite(ah))
    {
        double ahd = mag_get_d(ah);

        gi = sdq_choose(&X, &nreq, ctx, pc->pts[0], ahd, ahd, pc->top, logtol);
        if (gi >= 0)
        {
            /* the chord in the frame */
            acb_set(a, ctx->start);
            if (pc->fr->rel)
                acb_sub(a, a, &pc->fr->base, ctx->hp);
            sdq_get_cc(cc, ctx, pc->c);
            acb_set_cd(res, pc->pts[0]);
            sdq_chord_acb(res, ctx, cc, pc->top, a, res, pc->pts[0], pc->pts[0], gi);
            _sd_piece_c(cc, pc);
            sdq_gl_error(err, ctx, cc, b, ah, ah, X, sdq_steps[gi]);
            acb_add_error_mag(res, err);
            ok = 1;
        }
    }

    acb_clear(a);
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
    double * rephi, * mb, X, X2, nreq, nreq2, g[2], top = pc->top, sum;
    slong bprec = P->bprec;
    mag_t err, ah, zero;
    acb_t t, c, m;

    acb_zero(res);
    mag_init(err);
    mag_init(ah);
    mag_init(zero);
    acb_init(t);
    acb_init(c);
    acb_init(m);

    _sd_piece_c(c, pc);

    rephi = flint_malloc(sizeof(double) * 2 * (npts + 1));
    mb = rephi + npts + 1;

    /* estimates of the integrand at the vertices and of the integrals
       over the panels */
    for (i = 0; i < npts; i++)
    {
        cd phi, d1, d2;
        _sd_phi(&phi, &d1, &d2, P, pc->fr, pts[i]);
        rephi[i] = phi.re - top;
        if (i + 1 < npts)
            mb[i] = _sd_panel(g, P, pc->fr, pts[i], pts[i + 1], top);
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
            _sd_piece_vertex(m, pc, i);
            _sd_piece_vertex(t, pc, j);
            _sd_bound_segment(err, 1, ctx->alpha, ctx->sym, ctx->nu, m, t, c, bprec);
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

        sdq_chord(t, ctx, pc, i, j, gi, lmax);

        /* the error of the rule: the midpoint and the half-length of
           the chord between two vertices */
        _sd_piece_vertex(m, pc, i);
        {
            acb_t b;
            acb_init(b);
            _sd_piece_vertex(b, pc, j);
            acb_sub(b, b, m, bprec);
            acb_mul_2exp_si(b, b, -1);
            acb_get_mag(ah, b);
            acb_add(m, m, b, bprec);
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
    /* (from the exponent: mag_get_d_log2_approx is clamped to the
       range of a small fmpz, 2^30 on 32-bit systems) */
    if (mag_is_zero(t))
        r = (double) COEFF_MIN * SD_LN2;
    else if (mag_is_inf(t))
        r = (double) COEFF_MAX * SD_LN2;
    else
        r = (fmpz_get_d(MAG_EXPREF(t)) + log2(MAG_MAN(t) * ldexp(1.0, -MAG_BITS))) * SD_LN2;
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
    double logL = C->logL[0], delta, cutoff, mag, wmax, xlo, xhi, qprec, relmag = 0.0;
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
            _sd_phi(&phi, &d1, &d2, P, pc->fr, pc->pts[k]);
            if (k == 0 || phi.re - pc->top > -cutoff)
            {
                wmax = fmax(wmax, cd_abs(pc->pts[k]));
                xlo = fmin(xlo, pc->pts[k].re);
                xhi = fmax(xhi, pc->pts[k].re);
                /* (the terms of the exponent in a frame with rel = 1) */
                if (pc->fr->rel)
                    relmag = fmax(relmag, _sd_relsize(pc->fr, pc->pts[k]) + fabs(pc->top));
            }
        }
    }

    /* the cancellation in phi(w) - c */
    mag = P->aal * exp(xhi) + ((P->abe != 0.0) ? P->abe * exp(-xlo) : 0.0)
        + P->anu * wmax + cd_abs(cref) + 2.0;
    if (C->big)
        mag = relmag + 2.0;
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
    /* (In a frame with rel = 1, the vertices have no rounding errors
       relative to the width of the saddle, and the errors of the
       terms of the exponent are relative to their size.) */
    if (C->big)
        ctx.bits = qprec + SDQ_GUARD + log2(mag);
    ctx.cabs = cd_abs(cref);
    ctx.wp = (slong) ctx.bits + 8;
    ctx.F = NULL;
    ctx.kF = NULL;
    ctx.hp = ctx.wp + P->bprec;
    if (C->big)
    {
        acb_init(&ctx.ka); acb_init(&ctx.kb); acb_init(&ctx.kr); acb_init(&ctx.kk);
        acb_init(&ctx.ks); acb_init(&ctx.kd);
        acb_init(&ctx.shift);
        ctx.coef = _acb_vec_init(SDQ_REL_NC);
    }
    ctx.par_set[0] = ctx.par_set[1] = ctx.par_set[2] = ctx.par_set[3] = 0;
    ctx.budget = 4000 + 40 * (slong) qprec;
    ctx.scratch = NULL;

    if (dfloat_is_supported() && ctx.bits <= 4 * 53)
    {
        ctx.N = (int) ceil(fmax(ctx.bits / 53.0, 1.0));
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
        sdq_set_frame(&ctx, pc);
        if (C->big)
            ctx.cabs = fabs(pc->top);

        ok = sdq_piece_integrate(t, &ctx, pc, logtol, cutoff);

        if (ok && i == 0 && start != NULL)
        {
            ok = sdq_start_chord(c, &ctx, pc, logtol);
            /* (the piece is traversed from the starting point: sign = 1) */
            acb_add(t, t, c, ctx.wp);
        }

        if (pc->sign < 0)
            acb_neg(t, t);

        if (C->big)
        {
            acb_sub(c, pc->cx, &C->cref, ctx.hp);
        }
        else
        {
            arb_set_d(acb_realref(c), pc->c.re);
            arf_set_d(u, cref.re);
            arb_sub_arf(acb_realref(c), acb_realref(c), u, ctx.wp);
            arb_set_d(acb_imagref(c), pc->c.im);
            arf_set_d(u, cref.im);
            arb_sub_arf(acb_imagref(c), acb_imagref(c), u, ctx.wp);
        }
        /* (acb_exp gives up on numbers beyond 2^(2 prec), and the
           bounds of a negligible piece can be as large) */
        acb_exp(c, c, C->big ? FLINT_MAX(ctx.wp, P->bprec) : ctx.wp);
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
        if (C->big)
            acb_set(c, &C->cref);
        else
            acb_set_cd(c, cref);
        /* (acb_exp gives up on numbers beyond 2^(2 prec)) */
        acb_exp(c, c, C->big ? FLINT_MAX(ctx.wp, P->bprec) : ctx.wp);
        acb_mul(res, sum, c, ctx.wp);

        ok = acb_is_finite(res);
    }

    flint_free(ctx.scratch);
    if (C->big)
    {
        acb_clear(&ctx.ka); acb_clear(&ctx.kb); acb_clear(&ctx.kr); acb_clear(&ctx.kk);
        acb_clear(&ctx.ks); acb_clear(&ctx.kd);
        acb_clear(&ctx.shift);
        _acb_vec_clear(ctx.coef, SDQ_REL_NC);
    }
    acb_clear(sum);
    acb_clear(gsum);
    acb_clear(t);
    acb_clear(c);
    arf_clear(u);
    return ok;
}
