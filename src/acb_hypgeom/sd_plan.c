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
    Steepest-descent contours for integrals of exp(phi(w)),

        phi(w) = alpha e^w + beta e^(-w) - nu w,

    between two valleys at infinity, or from a finite starting point to
    a valley.

    EVERYTHING IN THIS FILE IS HEURISTIC.  It is done in double
    precision with the functions of the C library, without regard to
    rounding errors, overflow or underflow, and nothing computed here
    enters an error bound.  Since the integrand is entire, any path
    between the two ends gives the same integral, so the contour found
    here only determines how accurate and how expensive the evaluation
    along it is.  The output consists of exact data: polygons with
    vertices in double precision, reference values c in double
    precision, and the kinds of the ends of the polygons (a valley, with
    its index, or none).  The rigorous part is in sd_bound.c: bounds
    for the integrand on disks, for its integral over segments, and for
    the tails along rays from the ends of the polygons into the valleys
    (where it is also verified that a ray stays in the valley claimed
    here), all for parameters in the input balls, in arb and mag
    arithmetic.  The quantities called bounds below (the estimates of
    the L1 norms logL, logtail and the panel estimates) are the double
    precision counterparts of those bounds; they are used to decide
    where a path can end, to choose between routes, and to set the
    tolerances of the quadrature and of the expansions.

    Along a steepest-descent path the integrand has constant phase, so
    the L1 norm is close to the absolute value of the integral for a
    single saddle, and the integral can be evaluated along the path
    without cancellation (sd_integrate.c, sd_expand.c).

    The contour is planned as follows: the saddle points of phi (the
    zeros of phi', given by a quadratic equation in e^w, with a pair of
    nearby saddles replaced by a single cubic one), the descent paths
    from them traced by the midpoint rule (ending in a valley or, on a
    Stokes line, in a lower saddle), and a search for the cheapest
    route through the resulting graph.  A finite starting point is
    treated as a saddle with a single descent path.  On each straight
    panel, Re phi is majorized by its second-order Taylor polynomial at
    the panel midpoint plus |h|^3 sup |phi'''| / 6 (|phi'''| <=
    |alpha| e^x + |beta| e^(-x)), whose exponential integrates in closed
    form; the tails from the ends of the polygon into the valleys are
    estimated along rays in closed form.  Where the descent flow does
    not separate the valleys (flat landscapes: small parameters), a
    fixed polygon of straight segments is used instead.  The parameters
    must be below 2^44 or so.
*/

#include "acb.h"
#include "acb_hypgeom/sd.h"

static cd cd_div(cd a, cd b)
{
    double d = b.re * b.re + b.im * b.im;
    return cd_make((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d);
}

static cd cd_sqrt(cd a)
{
    double r = cd_abs(a), re, im;
    re = sqrt(0.5 * (r + a.re));
    im = sqrt(0.5 * fmax(r - a.re, 0.0));
    if (a.im < 0) im = -im;
    return cd_make(re, im);
}

static cd cd_log(cd a) { return cd_make(log(cd_abs(a)), cd_arg(a)); }

/* phi, phi', phi'' at w; phi''' = phi' + nu */
void
_sd_phi(cd * phi, cd * d1, cd * d2, const sd_problem * P, cd w)
{
    double e = exp(w.re), sn = sin(w.im), cs = cos(w.im);
    cd a = cd_mul(P->al, cd_make(e * cs, e * sn)), b = cd_make(0.0, 0.0);

    if (P->abe != 0.0)
        b = cd_mul(P->be, cd_make(cs / e, -sn / e));

    *d2 = cd_add(a, b);
    *phi = cd_sub(*d2, cd_mul(P->nu, w));
    *d1 = cd_sub(cd_sub(a, b), P->nu);
}

/* |alpha| e^xhi + |beta| e^(-xlo) */
static double
sd_size(const sd_problem * P, double xlo, double xhi)
{
    return P->aal * exp(xhi) + ((P->abe != 0.0) ? P->abe * exp(-xlo) : 0.0);
}

/* r_alpha e^xhi + r_beta e^(-xlo) */
static double
sd_ball(const sd_problem * P, double xlo, double xhi)
{
    return P->ral * exp(xhi) + ((P->abe != 0.0) ? P->rbe * exp(-xlo) : 0.0);
}

/* the weight e^xhi + [beta != 0] e^(-xlo) */
static double
sd_weight(const sd_problem * P, double xlo, double xhi)
{
    return exp(xhi) + ((P->abe != 0.0) ? exp(-xlo) : 0.0);
}

/* an allowance for the rounding of Re phi evaluated in double
   precision at w (in the estimates) */
double
_sd_margin(const sd_problem * P, cd w)
{
    return 16.0 * SD_EPS * (sd_size(P, w.re, w.re) + P->anu * cd_abs(w) + 1.0);
}

double
_sd_logaddexp(double a, double b)
{
    if (a == -INFINITY) return b;
    if (b == -INFINITY) return a;
    return fmax(a, b) + log1p(exp(-fabs(a - b)));
}

static int
sd_end_eq(sd_end a, sd_end b)
{
    return a.type == b.type && (a.type == SD_LINEAR || a.type == SD_START || a.k == b.k);
}

/* the centre of a valley (the imaginary part) */
static double
sd_centre(const sd_problem * P, sd_end v)
{
    if (v.type == SD_PLUS) return SD_PI - P->arg_al + 2 * SD_PI * v.k;
    if (v.type == SD_MINUS) return P->arg_mbe + 2 * SD_PI * v.k;
    if (v.type == SD_START) return P->start.im;
    return 20.0 * P->lin_dir.im;
}

/* log of int_1^inf exp(-K u) u^a du (an estimate from above), K > 0 */
static double
sd_log_gamma_tail(double K, double a)
{
    double v;

    if (a <= 0.0)
        return -K - log(K);
    /* u^a e^{-Ku/2} <= max(e^{-K/2}, (2a/(eK))^a), times int e^{-Ku/2} */
    v = log(2.0 / K) - 0.5 * K + fmax(-0.5 * K, a * log(2.0 * a / (2.718281828459045 * K)));
    /* u^a <= e^{a(u-1)} */
    if (a < K)
        v = fmin(v, -K - log(K - a));
    return v;
}

/* Estimates of the tail bounds (logarithms, relative to top; the
   rigorous counterpart is _sd_bound_tail) along the
   horizontal ray from we to Re w -> +inf, for the three weights;
   returns 0 if the ray is not certifiably downhill.  With w = x + i ye
   and u = e^(x - xe) >= 1,

       Re(alpha e^w) = -|alpha| c e^xe u,    c = -cos(arg alpha + ye),
       Re(beta e^-w) = Re(beta e^-we) / u,

   so with the ball terms the exponent is at most G - K (u - 1) + a log u,
   where G is its value at we, K = (|alpha| c - r_alpha) e^xe and
   a = max(-Re(beta e^-we), 0) + max(-Re nu, 0) + r_nu, using
   1 - 1/u <= log u and x - xe = log u. */
static int
sd_tail_plus(double tail[3], const sd_problem * P, cd we, double top, int * k)
{
    double c, K, a, G, xe = we.re, ye = we.im, e, ex = exp(xe), bt = 0.0;
    cd phi, d1, d2;

    c = -cos(P->arg_al + ye);
    if (!(c > 0.0))
        return 0;

    K = (P->aal * c - P->ral) * ex;
    if (!(K > 0.0))
        return 0;

    if (P->abe != 0.0)
        bt = cd_mul(P->be, cd_make(cos(ye) / ex, -sin(ye) / ex)).re;

    a = fmax(-bt, 0.0) + fmax(-P->nu.re, 0.0) + P->rnu;

    _sd_phi(&phi, &d1, &d2, P, we);
    G = phi.re + sd_ball(P, xe, xe) + P->rnu * cd_abs(we) + 2 * _sd_margin(P, we);

    /* int_xe^inf ... dx = int_1^inf e^G e^(-K(u-1)) u^a du / u */
    e = G - top + K;
    tail[0] = e + sd_log_gamma_tail(K, a - 1.0);
    tail[1] = e + log(sd_weight(P, xe, xe)) + sd_log_gamma_tail(K, a);
    tail[2] = e + log(fabs(xe) + fabs(ye) + 1.0) + sd_log_gamma_tail(K, a);
    *k = (int) floor((ye + P->arg_al - SD_PI) / (2 * SD_PI) + 0.5);
    return 1;
}

/* the ray to Re w -> -inf: phi(-w; alpha, beta, nu) = phi(w; beta, alpha, -nu) */
static int
sd_tail_minus(double tail[3], const sd_problem * P, cd we, double top, int * k)
{
    sd_problem Q = *P;
    int ok;

    if (P->abe == 0.0)
        return 0;

    Q.al = P->be; Q.be = P->al;
    Q.aal = P->abe; Q.abe = P->aal;
    Q.ral = P->rbe; Q.rbe = P->ral;
    Q.nu = cd_neg(P->nu);
    Q.arg_al = cd_arg(Q.al);
    ok = sd_tail_plus(tail, &Q, cd_neg(we), top, k);
    if (ok)
        *k = (int) floor((we.im - P->arg_mbe) / (2 * SD_PI) + 0.5);
    return ok;
}

/* the ray we + t e^(i theta), t >= 0, in the direction of the linear
   valley: Re(nu w) grows at the rate mu > 0, and with cos(theta) <= 0
   (cos(theta) = 0 if beta != 0), |alpha e^w + beta e^-w| stays below
   its bound at we */
static int
sd_tail_linear(double tail[3], const sd_problem * P, cd we, double top, int * k)
{
    double mu = P->lin_mu, G, e;

    if (!(mu > 0.0))
        return 0;

    G = sd_size(P, we.re, we.re) + sd_ball(P, we.re, we.re)
        - cd_mul(P->nu, we).re + P->rnu * cd_abs(we) + 2 * _sd_margin(P, we);

    e = G - top;
    tail[0] = e - log(mu);
    tail[1] = e + log(sd_weight(P, we.re, we.re) / mu);
    tail[2] = e + log(cd_abs(we) / mu + 1.0 / (mu * mu));
    *k = 0;
    return 1;
}

static int
sd_tail(double tail[3], const sd_problem * P, cd we, double top, int type, int * k)
{
    if (type == SD_PLUS) return sd_tail_plus(tail, P, we, top, k);
    if (type == SD_MINUS) return sd_tail_minus(tail, P, we, top, k);
    return sd_tail_linear(tail, P, we, top, k);
}

/* The panel [wa, wb]: an estimate of the logarithm of the bound of
   _sd_bound_segment for the integral of exp(Re phi - top) |dw| over
   the panel and the input balls; g is set to the suprema of the second
   and third weights on the panel. */
double
_sd_panel(double g[2], const sd_problem * P, cd wa, cd wb, double top)
{
    cd h = cd_scale(cd_sub(wb, wa), 0.5), wm = cd_scale(cd_add(wa, wb), 0.5);
    cd phi, d1, d2, t;
    double A, B, C, xlo, xhi, wmax, ah, m, lv;

    _sd_phi(&phi, &d1, &d2, P, wm);

    ah = cd_abs(h);
    xlo = fmin(wa.re, wb.re);
    xhi = fmax(wa.re, wb.re);
    wmax = fmax(cd_abs(wa), cd_abs(wb));

    t = cd_mul(d1, h);
    B = t.re;
    t = cd_mul(d2, cd_mul(h, h));
    C = 0.5 * t.re;
    A = phi.re - top + ah * ah * ah * sd_size(P, xlo, xhi) / 6.0
        + sd_ball(P, xlo, xhi) + P->rnu * wmax
        + _sd_margin(P, wm) + 16.0 * SD_EPS * (fabs(B) + fabs(C));

    /* int_{-1}^{1} exp(A + B t + C t^2) dt */
    m = fmax(A - B + C, A + B + C);
    if (C < 0 && fabs(B) < -2 * C)
        m = fmax(m, A - B * B / (4 * C));
    lv = log(2.0) + m;
    if (C < 0)
        lv = fmin(lv, 0.5 * log(SD_PI / (-C)) + A - B * B / (4 * C));
    lv += log(ah) + 1e-12;

    g[0] = sd_weight(P, xlo, xhi);
    g[1] = wmax;
    return lv;
}

/* add the panel bound to the three sums */
static void
sd_panel_acc(double acc[3], const sd_problem * P, cd wa, cd wb, double top)
{
    double g[2], v;

    v = exp(_sd_panel(g, P, wa, wb, top));
    /* an upper bound also after underflow */
    if (!(v >= 1e-300))
        v = 1e-300;
    acc[0] += v;
    acc[1] += v * g[0];
    acc[2] += v * g[1];
}

/* the step length for following a path at w: resolve the local width
   and keep the cubic remainder of the panel majorant small */
static double
sd_step(const sd_problem * P, double xlo, double xhi, cd d2)
{
    double h = fmin(1.0, 0.5 / sqrt(cd_abs(d2) + 1e-300));
    return fmin(h, 1.0 / cbrt(sd_size(P, xlo, xhi) + 1e-300));
}

static void
sd_half_push(sd_half * H, cd w)
{
    if (H->npts < SD_MAXSTEPS + 2)
        H->pts[H->npts++] = w;
}

/* Trace the descent path from the saddle self in the direction dir,
   bounding the L1 norms along the polygon and the tail.  The path ends
   when a tail bound certifies a valley, or in a lower saddle; on
   failure, the label has type 0. */
static void
sd_trace(sd_half * H, const sd_problem * P, const sd_work * W, int self, cd dir)
{
    const sd_saddle * S = W->saddles;
    cd w = S[self].w, phi, d1, d2, d = dir, wm, wn;
    double acc[3] = { 0.0, 0.0, 0.0 }, top = S[self].c.re, h, S3, drop;
    int ntest = -100, step, i, j;

    H->saddle = self;
    H->label = sd_end_make(0, 0);
    H->pts = flint_malloc(sizeof(cd) * (SD_MAXSTEPS + 2));
    H->npts = 0;
    sd_half_push(H, w);

    _sd_phi(&phi, &d1, &d2, P, w);

    for (step = 0; step < SD_MAXSTEPS; step++)
    {
        /* give up outside the range of double precision, or when the
           descent is too slow to get anywhere */
        if (fabs(w.re) > SD_XMAX || fabs(w.im) > 1e6 || !(acc[1] < 1e280)
            || (step >= 50 && (top - phi.re) * SD_MAXSTEPS < P->drop * step))
            return;

        if (step != 0 && cd_abs(d1) > 0)
            d = cd_scale(cd_conj(d1), -1.0 / cd_abs(d1));

        /* The step length: resolve the local width, limit the drop per
           step and keep the cubic remainder of the panel majorant
           small.  Once the integrand is negligible, below e^-drop
           times the L1 norm so far (the polygon only continues until a
           tail can be certified), the panel bounds may be loose by
           e^slack, with slack growing with the drop. */
        S3 = sd_size(P, w.re, w.re) + 1e-300;
        drop = top - phi.re - (P->drop - log(acc[0]) + 4.0);
        if (drop < 0.0)
        {
            h = fmin(1.0, 1.0 / sqrt(cd_abs(d2) + 1e-300));
            if (cd_abs(d1) > 0)
                h = fmin(h, P->maxdrop / cd_abs(d1));
            h = fmin(h, 1.0 / cbrt(S3));
        }
        else
        {
            double slack = 2.0 + 0.5 * drop;
            h = fmin(1.0, sqrt(2 * slack / (cd_abs(d2) + 1e-300)));
            h = fmin(h, cbrt(6 * slack / S3));
        }
        if (step == 0)
            h = fmin(h, S[self].width);

        /* do not step over a lower saddle */
        for (i = 0; i < W->nsaddles; i++)
        {
            if (i != self && S[i].kind != 1 && S[i].c.re < top)
            {
                double dist = cd_abs(cd_sub(w, S[i].w));
                if (dist < 2.0)
                    h = fmin(h, fmax(0.5 * dist, 0.5 * S[i].width));
            }
        }

        /* midpoint rule for the direction (but the first step from a
           degenerate saddle is in the given direction) */
        wm = cd_add(w, cd_scale(d, 0.5 * h));
        _sd_phi(&phi, &d1, &d2, P, wm);
        if (cd_abs(d1) > 0 && !(step == 0 && S[self].kind == 3))
            d = cd_scale(cd_conj(d1), -1.0 / cd_abs(d1));
        wn = cd_add(w, cd_scale(d, h));

        sd_panel_acc(acc, P, w, wn, top);
        w = wn;
        sd_half_push(H, w);

        _sd_phi(&phi, &d1, &d2, P, w);
        drop = top - phi.re;

        /* running into a lower saddle: end there (the route continues
           along that saddle's own descent paths) */
        if (step >= 2)
        {
            for (i = 0; i < W->nsaddles; i++)
            {
                if (i != self && S[i].kind != 1
                    && cd_abs(cd_sub(w, S[i].w)) < 2.0 * S[i].width
                    && S[i].c.re < top - 1e-9 * (fabs(top) + 1.0))
                {
                    sd_panel_acc(acc, P, w, S[i].w, top);
                    if (!(acc[1] < 1e280))
                        return;
                    sd_half_push(H, S[i].w);
                    H->label = sd_end_make(SD_SADDLE, i);
                    for (j = 0; j < 3; j++)
                    {
                        H->logL[j] = log(acc[j]);
                        H->logtail[j] = -INFINITY;
                    }
                    return;
                }
            }
        }

        /* Deep in a valley, the tails along the rays of several
           valleys can be negligible; the path is heading into the one
           whose ray is best aligned with it. */
        if (drop > P->drop - log(acc[0]) && (step - ntest) >= 1)
        {
            double tail[3], best = -INFINITY, dot;
            int k;

            ntest = step;
            if (cd_abs(d1) > 0)
                d = cd_scale(cd_conj(d1), -1.0 / cd_abs(d1));

            for (i = SD_PLUS; i <= SD_LINEAR; i++)
            {
                dot = (i == SD_PLUS) ? d.re : (i == SD_MINUS) ? -d.re
                    : d.re * P->lin_dir.re + d.im * P->lin_dir.im;

                if (dot > best && sd_tail(tail, P, w, top, i, &k)
                    && tail[0] < log(acc[0]) - P->tailtol
                    && tail[1] < log(acc[1]) - P->tailtol
                    && tail[2] < log(acc[2]) - P->tailtol)
                {
                    best = dot;
                    H->label = sd_end_make(i, k);
                    for (j = 0; j < 3; j++)
                    {
                        H->logL[j] = _sd_logaddexp(log(acc[j]), tail[j]);
                        H->logtail[j] = tail[j];
                    }
                }
            }

            if (H->label.type != 0)
                return;
        }
    }
}

/* trace the descent paths of a saddle */
static void
sd_work_trace(sd_work * W, const sd_problem * P, int idx)
{
    sd_saddle * S = W->saddles + idx;
    cd dirs[3], phi, d1, d2;
    int m, ndir = S->kind;

    _sd_phi(&phi, &d1, &d2, P, S->w);

    if (S->kind == 1 && cd_abs(d1) > 1e-8 * (cd_abs(d2) + 1e-300))
    {
        dirs[0] = cd_scale(cd_conj(d1), -1.0 / cd_abs(d1));
    }
    else if (S->kind == 3)
    {
        /* phi''' = phi' + nu */
        double th = (SD_PI - cd_arg(cd_add(d1, P->nu))) / 3.0;
        for (m = 0; m < 3; m++)
            dirs[m] = cd_expi(th + 2 * SD_PI * m / 3.0);
    }
    else
    {
        dirs[0] = cd_expi((SD_PI - cd_arg(d2)) / 2.0);
        dirs[1] = cd_neg(dirs[0]);
    }

    for (m = 0; m < ndir && W->nhalf < SD_MAXHALF; m++)
    {
        sd_trace(W->halves + W->nhalf, P, W, idx, dirs[m]);
        S->half[S->nhalf++] = W->nhalf;
        W->nhalf++;
    }
    S->traced = 1;
}

static void
sd_work_add(sd_work * W, const sd_problem * P, cd ws, int kind, double sep)
{
    sd_saddle * S;
    cd d1, d2;

    if (W->nsaddles >= SD_MAXSADDLES)
        return;

    S = W->saddles + W->nsaddles++;
    S->w = ws;
    S->kind = kind;
    S->nhalf = 0;
    S->traced = 0;
    _sd_phi(&S->c, &d1, &d2, P, ws);
    /* the width of the saddle; all the derivatives of phi from the
       second are of the same order, so this is at most 1 */
    S->width = 1.0;
    if (kind == 3)
        S->width = fmin(1.0, fmax(0.5 / cbrt(cd_abs(cd_add(d1, P->nu)) + 1e-300), 3.0 * sep));
    else if (kind == 2)
        S->width = fmin(1.0, 0.5 / sqrt(cd_abs(d2) + 1e-300));
}

/* the saddles with |Im w| <= SD_YMAX are considered (the ends of the
   contours of interest are within 4 pi of the real axis) */
#define SD_YMAX (5.5 * SD_PI)

/* The saddles: with u = e^w, phi' = 0 is alpha u^2 - nu u - beta = 0.
   A pair of roots within the width |phi'''|^(-1/3) (phi''' = nu at a
   double root) of the resulting cubic saddle is replaced by that
   degenerate saddle at u = nu / (2 alpha) (if the roots are close).
   This point is not exactly a saddle: the descent paths from it start
   with a step along the three descent directions of the cubic term,
   beyond the two saddles, where that term determines the descent. */
void
_sd_work_init(sd_work * W, const sd_problem * P)
{
    cd u[2], base;
    double sep = 0.0;
    int nroots = 0, kind = 2, j, k;

    W->nsaddles = W->nhalf = 0;

    if (P->has_start)
        sd_work_add(W, P, P->start, 1, 0.0);

    if (P->abe == 0.0)
    {
        if (P->anu > 0.0)
            u[nroots++] = cd_div(P->nu, P->al);
    }
    else
    {
        cd sq = cd_sqrt(cd_add(cd_mul(P->nu, P->nu), cd_scale(cd_mul(P->al, P->be), 4.0))), t;

        /* nu + sq without cancellation; the product of the roots is -beta / alpha */
        if (cd_mul(P->nu, cd_conj(sq)).re < 0)
            sq = cd_neg(sq);
        t = cd_add(P->nu, sq);
        u[0] = cd_div(t, cd_scale(P->al, 2.0));
        u[1] = cd_div(cd_scale(P->be, -2.0), t);
        nroots = 2;

        sep = 0.5 * cd_abs(cd_log(cd_div(u[0], u[1])));
        if (sep < 0.5 * fmin(1.0, 1.0 / cbrt(P->anu)))
        {
            u[0] = cd_div(P->nu, cd_scale(P->al, 2.0));
            nroots = 1;
            kind = 3;
        }
    }

    for (j = 0; j < nroots; j++)
    {
        base = cd_log(u[j]);
        for (k = -3; k <= 3; k++)
        {
            cd ws = cd_make(base.re, base.im + 2 * SD_PI * k);
            if (fabs(ws.im) <= SD_YMAX && fabs(ws.re) <= SD_XMAX)
                sd_work_add(W, P, ws, kind, sep);
        }
    }
}

void
_sd_work_clear(sd_work * W)
{
    int i;
    for (i = 0; i < W->nhalf; i++)
        flint_free(W->halves[i].pts);
}

/* A route between two ends: a sequence of half-paths joined at saddles
   (each traversed in either direction), or a fixed polygon */
typedef struct
{
    double logL[3];          /* log of the L1 bounds (absolute) */
    int nedges;
    int half[SD_MAXROUTE];
    int sign[SD_MAXROUTE];   /* +1: traversed from its saddle outward */
    cd * pts;                /* the fixed polygon, with a tail at each end */
    int npts;
    double top;
    double logtail[2][3];
    sd_end ends[2];          /* of the fixed polygon */
}
sd_route;

/* Route search.  The saddles and the valleys are the nodes of a graph
   whose edges are the half-paths, with the L1 bounds as lengths: any
   path in the graph is a contour, and the shortest one from vs to ve is
   found by the Bellman-Ford algorithm. */
#define SD_MAXNODES (SD_MAXSADDLES + SD_MAXHALF + 2)

static int
sd_node(sd_end * nodes, int * n, sd_end v)
{
    int i;

    if (v.type == SD_SADDLE)
        return v.k;
    if (v.type == SD_START)
        return 0;
    for (i = SD_MAXSADDLES; i < *n; i++)
        if (sd_end_eq(nodes[i], v))
            return i;
    nodes[*n] = v;
    return (*n)++;
}

static int
sd_search(sd_route * R, const sd_work * W, sd_end vs, sd_end ve)
{
    sd_end nodes[SD_MAXNODES];
    double dist[SD_MAXNODES][3], d[3];
    int from[SD_MAXNODES], half[SD_MAXNODES], sign[SD_MAXNODES];
    int n = SD_MAXSADDLES, src, dst, round, changed, hi, k, u, v, dir, len;

    src = sd_node(nodes, &n, vs);
    dst = sd_node(nodes, &n, ve);
    for (hi = 0; hi < W->nhalf; hi++)
        if (W->halves[hi].label.type != 0)
            sd_node(nodes, &n, W->halves[hi].label);

    for (u = 0; u < n; u++)
    {
        dist[u][0] = INFINITY;
        from[u] = -1;
    }
    dist[src][0] = dist[src][1] = dist[src][2] = -INFINITY;

    for (round = 0, changed = 1; round < n && changed; round++)
    {
        changed = 0;

        for (hi = 0; hi < W->nhalf; hi++)
        {
            const sd_half * h = W->halves + hi;

            if (h->label.type == 0)
                continue;

            /* dir = 1: from the saddle of the half-path outward */
            for (dir = 1; dir >= -1; dir -= 2)
            {
                u = h->saddle;
                v = sd_node(nodes, &n, h->label);
                if (dir == -1)
                {
                    k = u; u = v; v = k;
                }

                /* the starting point is only left */
                if (dist[u][0] == INFINITY || (W->saddles[0].kind == 1 && v == 0))
                    continue;

                for (k = 0; k < 3; k++)
                    d[k] = _sd_logaddexp(dist[u][k], W->saddles[h->saddle].c.re + h->logL[k]);

                if (d[0] < dist[v][0] - 1e-9)
                {
                    for (k = 0; k < 3; k++)
                        dist[v][k] = d[k];
                    from[v] = u;
                    half[v] = hi;
                    sign[v] = dir;
                    changed = 1;
                }
            }
        }
    }

    if (dist[dst][0] == INFINITY)
        return 0;

    for (len = 0, v = dst; v != src; v = from[v])
        if (++len > SD_MAXROUTE)
            return 0;

    R->nedges = len;
    R->pts = NULL;
    for (k = 0; k < 3; k++)
        R->logL[k] = dist[dst][k];
    for (v = dst; v != src; v = from[v])
    {
        len--;
        R->half[len] = half[v];
        R->sign[len] = sign[v];
    }
    return 1;
}

/* The cheapest steepest-descent route from vs to ve; returns 0 if no
   route is found.  The saddles are traced in order of increasing
   Re phi (and closeness to the ends), stopping as soon as a route is
   found which no untraced saddle can improve. */
static int
sd_route_sd(sd_route * R, sd_work * W, const sd_problem * P, sd_end vs, sd_end ve)
{
    int k, i, found, order[SD_MAXSADDLES], n;
    double ya = sd_centre(P, vs), yb = sd_centre(P, ve), ymid = 0.5 * (ya + yb);
    double key[SD_MAXSADDLES];

    for (i = 0; i < W->nsaddles; i++)
    {
        order[i] = i;
        key[i] = (W->saddles[i].kind == 1) ? -INFINITY
            : W->saddles[i].c.re + 1e-6 * fabs(W->saddles[i].w.im - ymid);
    }
    for (i = 1; i < W->nsaddles; i++)
        for (k = i; k > 0 && key[order[k]] < key[order[k - 1]]; k--)
        {
            int t = order[k]; order[k] = order[k - 1]; order[k - 1] = t;
        }

    found = sd_search(R, W, vs, ve);

    for (n = 0; n < W->nsaddles; n++)
    {
        int idx = order[n];

        if (found && R->logL[0] < W->saddles[idx].c.re - 20.0)
            break;
        if (W->saddles[idx].traced)
            continue;

        /* only the saddles between the ends, or near */
        if (W->saddles[idx].w.im < fmin(ya, yb) - 1.5 * SD_PI
            || W->saddles[idx].w.im > fmax(ya, yb) + 1.5 * SD_PI)
            continue;

        /* a saddle more than pi beyond those of the route found is a
           periodic copy of a nearer saddle (phi(w + 2 pi i) = phi(w)
           - 2 pi i nu): it cannot give a cheaper route */
        if (found)
        {
            double far = 0.0;
            for (i = 0; i < R->nedges; i++)
                far = fmax(far, fabs(W->saddles[W->halves[R->half[i]].saddle].w.im - ymid));
            if (fabs(W->saddles[idx].w.im - ymid) > far + SD_PI)
                continue;
        }

        sd_work_trace(W, P, idx);
        found = sd_search(R, W, vs, ve);
    }

    return found;
}

/* A fixed polygon through the vertices base[0..nb-1], where each end
   is the starting point or continues along a horizontal ray in the
   direction dir into the valley ends[i].  The rays are extended until
   their tails are certified and below e^-tailtol times the L1 norm of
   the polygon; the segments are subdivided into panels. */
#define SD_FIXED_MAXSTEPS 1000
static int
sd_fixed_build(sd_route * R, const sd_problem * P, const cd * base, int nb,
    const sd_end * ends, const double * dir)
{
    double ext[2] = { 0.0, 0.0 }, acc[3], tail[3], top = 0.0;
    cd p[6], * pts;
    int np, n = 0, i, k, attempt, extend[2];

    pts = flint_malloc(sizeof(cd) * (SD_FIXED_MAXSTEPS + 2));

    /* the tail of a horizontal ray is negligible about where
       |alpha| e^x (or |beta| e^-x) exceeds tailtol */
    for (i = 0; i < 2; i++)
    {
        const cd b = (i == 0) ? base[0] : base[nb - 1];
        if (ends[i].type == SD_PLUS)
            ext[i] = fmax(0.0, log(P->tailtol / P->aal) - b.re);
        if (ends[i].type == SD_MINUS)
            ext[i] = fmax(0.0, log(P->tailtol / P->abe) + b.re);
    }

    for (attempt = 0; attempt < 20; attempt++)
    {
        np = 0;
        if (ends[0].type != SD_START)
            p[np++] = cd_make(base[0].re + dir[0] * ext[0], base[0].im);
        for (i = 0; i < nb; i++)
            p[np++] = base[i];
        if (ends[1].type != SD_START)
            p[np++] = cd_make(base[nb - 1].re + dir[1] * ext[1], base[nb - 1].im);

        if (fabs(p[0].re) > SD_XMAX || fabs(p[np - 1].re) > SD_XMAX)
            break;

        /* subdivide */
        n = 0;
        pts[n++] = p[0];
        for (i = 0; i + 1 < np; i++)
        {
            cd a = p[i], b = p[i + 1], d = cd_sub(b, a), phi, d1, d2;
            double len = cd_abs(d), t = 0.0;

            while (t < len && n < SD_FIXED_MAXSTEPS)
            {
                _sd_phi(&phi, &d1, &d2, P, pts[n - 1]);
                t += sd_step(P, fmin(pts[n - 1].re, b.re), fmax(pts[n - 1].re, b.re), d2);
                pts[n++] = (t >= len) ? b : cd_add(a, cd_scale(d, t / len));
            }
            if (t < len)
                goto fail;
        }

        /* the reference level: the maximum of Re phi at the vertices */
        for (i = 0; i < n; i++)
        {
            cd phi, d1, d2;
            _sd_phi(&phi, &d1, &d2, P, pts[i]);
            top = (i == 0) ? phi.re : fmax(top, phi.re);
        }

        acc[0] = acc[1] = acc[2] = 0.0;
        for (i = 0; i + 1 < n; i++)
            sd_panel_acc(acc, P, pts[i], pts[i + 1], top);

        extend[0] = extend[1] = 0;
        for (i = 0; i < 2; i++)
        {
            for (k = 0; k < 3; k++)
                R->logtail[i][k] = -INFINITY;

            if (ends[i].type == SD_START)
                continue;

            if (!sd_tail(tail, P, (i == 0) ? pts[0] : pts[n - 1], top, ends[i].type, &k))
                extend[i] = 1;
            else if (!sd_end_eq(sd_end_make(ends[i].type, k), ends[i]))
                goto fail;
            else
            {
                for (k = 0; k < 3; k++)
                {
                    R->logtail[i][k] = tail[k];
                    if (!(acc[k] > 0.0 && tail[k] < log(acc[k]) - P->tailtol))
                        extend[i] = 1;
                }
            }

            if (extend[i])
                ext[i] += 1.0;
        }

        if (!extend[0] && !extend[1])
        {
            R->nedges = 0;
            R->pts = pts;
            R->npts = n;
            R->top = top;
            R->ends[0] = ends[0];
            R->ends[1] = ends[1];
            for (k = 0; k < 3; k++)
                R->logL[k] = top + _sd_logaddexp(log(acc[k]),
                    _sd_logaddexp(R->logtail[0][k], R->logtail[1][k]));
            return 1;
        }
    }

fail:
    flint_free(pts);
    return 0;
}

/* The classical polygons, as a fallback where the descent flow does
   not separate the valleys: between two valleys V+-(k) of the
   Bessel-type problem (beta != 0), rays along the valley centres
   starting where |alpha| e^x, |beta| e^(-x) >= 1, joined by straight
   segments; from the starting point of the gamma-type problem
   (beta = 0), vertically to the centre of the valley V+(k) and then
   along it. */
static int
sd_route_fixed(sd_route * R, const sd_problem * P, sd_end vs, sd_end ve)
{
    cd base[4];
    double dir[2] = { 0.0, 0.0 };
    sd_end ends[2];
    int nb = 0, i;

    ends[0] = vs;
    ends[1] = ve;

    if (vs.type == SD_START && ve.type == SD_PLUS && P->abe == 0.0)
    {
        double y = sd_centre(P, ve), x = P->start.re;

        base[nb++] = P->start;
        base[nb++] = cd_make(x, y);
        base[nb++] = cd_make(fmax(x, -log(P->aal)), y);
        dir[1] = 1.0;
    }
    else if ((vs.type == SD_PLUS || vs.type == SD_MINUS)
          && (ve.type == SD_PLUS || ve.type == SD_MINUS) && P->abe != 0.0)
    {
        double xr0 = fmax(-log(fmin(P->aal, P->abe)), 0.0);

        /* the base of the ray of each end, and for a ray to the left
           the point of the spine x = xr0 at its height */
        for (i = 0; i < 2; i++)
        {
            double y = sd_centre(P, ends[i]);

            dir[i] = (ends[i].type == SD_PLUS) ? 1.0 : -1.0;
            if (i == 1 && dir[i] < 0)
                base[nb++] = cd_make(xr0, y);
            base[nb++] = cd_make(dir[i] * xr0, y);
            if (i == 0 && dir[i] < 0)
                base[nb++] = cd_make(xr0, y);
        }
    }
    else
        return 0;

    return sd_fixed_build(R, P, base, nb, ends, dir);
}

static void
sd_route_clear(sd_route * R)
{
    if (R->nedges == 0)
        flint_free(R->pts);
}

/* add the pieces of a route to the contour, with the sign sgn, and
   divided by 1 - q if geom is set: with the upper bound qmax for |q|,
   |1 / (1 - q)| <= 1 / (1 - qmax), and its derivative with respect to nu
   is at most 2 pi qmax / (1 - qmax)^2 */
static void
sd_contour_add(sd_contour * C, const sd_route * R, const sd_work * W, int sgn, int geom, double qmax)
{
    sd_piece * pc;
    double f = geom ? -log1p(-qmax) : 0.0;
    int e, k, first = C->npieces;

    for (k = 0; k < 3; k++)
        C->logL[k] = _sd_logaddexp(C->logL[k], R->logL[k] + f);
    if (geom)
        C->logL[2] = _sd_logaddexp(C->logL[2], R->logL[0] + log(2 * SD_PI * qmax) + 2 * f);

    if (R->nedges == 0)
    {
        pc = C->piece + C->npieces++;
        pc->pts = R->pts;
        pc->npts = R->npts;
        pc->c = cd_make(R->top, 0.0);
        pc->sign = sgn;
        pc->logtail[0] = R->logtail[0][0];
        pc->logtail[1] = R->logtail[1][0];
        pc->end[0] = R->ends[0];
        pc->end[1] = R->ends[1];
        pc->mid = -1;
        C->own[C->nown++] = R->pts;
    }

    for (e = 0; e < R->nedges; e++)
    {
        const sd_half * h = W->halves + R->half[e], * h2 = h + 1;

        pc = C->piece + C->npieces++;
        pc->pts = h->pts;
        pc->npts = h->npts;
        pc->c = W->saddles[h->saddle].c;
        pc->sign = sgn * R->sign[e];
        pc->logtail[0] = -INFINITY;
        pc->logtail[1] = h->logtail[0];
        pc->end[0] = sd_end_make(0, 0);
        pc->end[1] = h->label;
        pc->mid = -1;

        /* passing through a saddle: join the two half-paths, so that
           the quadrature can use chords across the saddle */
        if (e + 1 < R->nedges && R->sign[e] == -1 && R->sign[e + 1] == 1
            && (h2 = W->halves + R->half[e + 1])->saddle == h->saddle)
        {
            cd * pts = flint_malloc(sizeof(cd) * (h->npts + h2->npts));

            for (k = 0; k < h->npts; k++)
                pts[k] = h->pts[h->npts - 1 - k];
            for (k = 1; k < h2->npts; k++)
                pts[h->npts - 1 + k] = h2->pts[k];

            C->own[C->nown++] = pts;
            pc->pts = pts;
            pc->npts = h->npts + h2->npts - 1;
            pc->sign = sgn;
            pc->logtail[0] = pc->logtail[1];
            pc->logtail[1] = h2->logtail[0];
            pc->end[0] = pc->end[1];
            pc->end[1] = h2->label;
            /* (for the saddle-point expansion) */
            if (W->saddles[h->saddle].kind == 2
                && h->label.type != SD_SADDLE && h2->label.type != SD_SADDLE)
                pc->mid = h->npts - 1;
            e++;
        }
    }

    for (e = first; e < C->npieces; e++)
        C->piece[e].geom = geom;
}

void
_sd_contour_clear(sd_contour * C)
{
    int i;
    for (i = 0; i < C->nown; i++)
        flint_free(C->own[i]);
}

/* The route from the valley v = V(k) (of either family) to the linear
   valley, in the vertical direction sigma i: a steepest-descent route
   if there is one.  Otherwise, the valleys V(k + j sigma), j -> inf,
   lead to the linear valley, separated from each other by periodic
   copies of the same saddles.  Since phi(w + 2 pi i sigma) = phi(w) -
   2 pi i sigma nu, the integral from V(k + sigma) to V(k + 2 sigma) is q
   times that from V(k) to V(k + sigma), q = exp(-2 pi i sigma nu),
   |q| < 1, so that

       int_{V(k)}^{linear valley} = (1 / (1 - q)) int_{V(k)}^{V(k + sigma)}.

   Returns 1 or 2 (with the factor) if a route is found. */
static int
sd_route_linear(sd_route * R, sd_work * W, const sd_problem * P, sd_end v)
{
    int sigma = (P->lin_dir.im > 0.0) ? 1 : -1;

    if (sd_route_sd(R, W, P, v, sd_end_make(SD_LINEAR, 0)))
        return 1;
    if (sd_route_sd(R, W, P, v, sd_end_make(v.type, v.k + sigma)))
        return 2;
    return 0;
}

/* The steepest-descent route from vs to ve if there is one; otherwise,
   for two valleys of the Bessel-type problem with a linear valley,
   the routes from each of them to the linear valley; failing that, or
   if its L1 bound is smaller, the fixed polygon if fixed is set.

   With a starting point, the L1 bounds include the segment from the
   true starting point to its approximation, bounded on the disk around
   the latter. */
int
_sd_plan(sd_contour * C, sd_work * W, const sd_problem * P, sd_end vs, sd_end ve, int fixed)
{
    sd_route R, A, B;
    int found, ga, gb, k;

    C->npieces = 0;
    C->nown = 0;
    C->logL[0] = C->logL[1] = C->logL[2] = -INFINITY;

    found = sd_route_sd(&R, W, P, vs, ve);

    if (!found)
    {
        if (fixed)
            found = sd_route_fixed(&R, P, vs, ve);

        if (P->abe != 0.0 && P->lin_mu > 0.0
            && (vs.type == SD_PLUS || vs.type == SD_MINUS)
            && (ve.type == SD_PLUS || ve.type == SD_MINUS)
            && (ga = sd_route_linear(&A, W, P, vs)) != 0
            && (gb = sd_route_linear(&B, W, P, ve)) != 0)
        {
            /* an upper bound for |q| */
            double q = exp(-2 * SD_PI * P->lin_mu);

            sd_contour_add(C, &A, W, 1, ga == 2, q);
            sd_contour_add(C, &B, W, -1, gb == 2, q);

            if (!found || C->logL[0] < R.logL[0])
            {
                if (found)
                    sd_route_clear(&R);
                return 1;
            }

            _sd_contour_clear(C);
            C->npieces = 0;
            C->nown = 0;
            C->logL[0] = C->logL[1] = C->logL[2] = -INFINITY;
        }
    }

    if (!found)
        return 0;

    sd_contour_add(C, &R, W, 1, 0, 0.0);

    if (vs.type == SD_START && P->start_rad > 0.0)
    {
        cd phi, d1, d2, w = P->start;
        double r = P->start_rad, sh = sinh(0.5 * r), L;

        /* (with cosh r - 1 and sinh r - r rounded up) */
        _sd_phi(&phi, &d1, &d2, P, w);
        L = log(r) + phi.re + cd_abs(d1) * r + cd_abs(d2) * 2.0 * sh * sh * (1.0 + 1e-12)
            + cd_abs(cd_add(d1, P->nu)) * ((r < 0.5) ? r * r * r * 0.17 : sinh(r) - r)
            + sd_ball(P, w.re - r, w.re + r) + P->rnu * (cd_abs(w) + r)
            + 32.0 * SD_EPS * (sd_size(P, w.re - r, w.re + r) + P->anu * (cd_abs(w) + r) + 1.0);

        for (k = 0; k < 3; k++)
            C->logL[k] = _sd_logaddexp(C->logL[k], L + ((k == 0) ? 0.0
                : (k == 1) ? log(sd_weight(P, w.re - r, w.re + r)) : log(cd_abs(w) + r)));
    }

    return 1;
}

/* res = exp(L) approximately (for tolerances and estimates only) */
void
_sd_mag_approx_exp_d(mag_t res, double L)
{
    double e;

    if (L != L || L > 1e15)
    {
        mag_inf(res);
        return;
    }
    if (L < -1e15)
    {
        mag_zero(res);
        return;
    }

    L += 1e-9 + 1e-12 * fabs(L);
    e = floor(L / SD_LN2);
    mag_set_d(res, exp(L - e * SD_LN2) * (1.0 + 1e-12));
    mag_mul_2exp_si(res, res, (slong) e);
}

/* the midpoint of x as a complex double and an upper bound for the
   distance from it to any point of x; returns 0 if out of range */
static int
sd_get_cd(cd * mid, double * rad, const acb_t x, slong maxexp)
{
    mag_t t;

    if (!acb_is_finite(x)
        || arf_cmpabs_2exp_si(arb_midref(acb_realref(x)), maxexp) >= 0
        || arf_cmpabs_2exp_si(arb_midref(acb_imagref(x)), maxexp) >= 0)
        return 0;

    mid->re = arf_get_d(arb_midref(acb_realref(x)), ARF_RND_NEAR);
    mid->im = arf_get_d(arb_midref(acb_imagref(x)), ARF_RND_NEAR);

    mag_init(t);
    mag_hypot(t, arb_radref(acb_realref(x)), arb_radref(acb_imagref(x)));
    *rad = mag_get_d(t) * (1.0 + 1e-12) + 4.0 * SD_EPS * cd_abs(*mid);
    mag_clear(t);
    return *rad < 1e300;
}

int
_sd_problem_init(sd_problem * P, const acb_t alpha, int sym, const acb_t nu,
    const acb_t start, double drop, double tailtol, double maxdrop)
{
    if (!sd_get_cd(&P->al, &P->ral, alpha, 44) || !sd_get_cd(&P->nu, &P->rnu, nu, 44))
        return 0;

    /* (for tiny alpha, the contours are long, and other methods apply) */
    P->aal = cd_abs(P->al);
    P->anu = cd_abs(P->nu);
    if (!(P->aal > 1e-120))
        return 0;

    if (sym)
    {
        P->be = cd_neg(P->al);
        P->rbe = P->ral;
        P->abe = P->aal;
    }
    else
    {
        P->be = cd_make(0.0, 0.0);
        P->rbe = P->abe = 0.0;
    }

    P->arg_al = cd_arg(P->al);
    P->arg_mbe = cd_arg(cd_neg(P->be));

    /* The linear valley: Re(nu e^(i theta)) is largest in the
       direction conj(nu), which is admissible if beta = 0 and it
       points to the left; otherwise the vertical direction. */
    if (!sym && P->nu.re <= 0.0 && P->anu > 0.0)
    {
        P->lin_dir = cd_scale(cd_conj(P->nu), 1.0 / P->anu);
        P->lin_mu = P->anu - P->rnu;
    }
    else
    {
        P->lin_dir = cd_make(0.0, (P->nu.im > 0.0) ? -1.0 : 1.0);
        P->lin_mu = fabs(P->nu.im) - P->rnu;
    }
    /* rounding in Re(nu e^(i theta)) */
    P->lin_mu -= 8.0 * SD_EPS * P->anu;

    P->has_start = (start != NULL);
    if (P->has_start && (!sd_get_cd(&P->start, &P->start_rad, start, 10)
        || fabs(P->start.re) > 280.0))
        return 0;

    P->drop = drop;
    P->tailtol = tailtol;
    P->maxdrop = maxdrop;
    return 1;
}
