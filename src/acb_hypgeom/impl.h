/*
    Copyright (C) 2025 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ACB_HYPGEOM_IMPL_H
#define ACB_HYPGEOM_IMPL_H

#include "acb_types.h"

int _mag_gt_norm_ui(const mag_t a, const mag_t b, const mag_t c, ulong n);

void acb_hypgeom_gamma_stirling_choose_param(int * reflect, slong * r, slong * n, const acb_t z, int use_reflect, int digamma, slong prec);
void acb_hypgeom_gamma_stirling_inner(acb_t s, const acb_t z, slong N, slong prec);
slong acb_hypgeom_pfq_choose_n_max(acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong prec, slong n_max);
slong _acb_hypgeom_pfq_choose_n_bound(mag_t C, int * have_bound, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong prec, slong n_max);
void _acb_hypgeom_pfq_direct_max(acb_t res, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong prec, slong n_max);
int acb_hypgeom_u_asymp_determine_region(const mag_t r, const mag_t zlo, const acb_t z);
void acb_hypgeom_mag_chi(mag_t chi, ulong n);
void _acb_hypgeom_legendre_q_single(acb_t res, const acb_t n, const acb_t m, const acb_t z, slong prec);
void _acb_hypgeom_legendre_q_double(acb_t res, const acb_t n, const acb_t m, const acb_t z, slong prec);

void acb_poly_reciprocal_majorant(arb_poly_t res, const acb_poly_t poly, slong prec);

/* Functions selected by an integer, for the Bessel integration code */
#define ACB_HYPGEOM_KIND_BESSEL_J 0
#define ACB_HYPGEOM_KIND_BESSEL_Y 1
#define ACB_HYPGEOM_KIND_BESSEL_JY 2
#define ACB_HYPGEOM_KIND_HANKEL_1 3
#define ACB_HYPGEOM_KIND_HANKEL_2 4
#define ACB_HYPGEOM_KIND_BESSEL_I 5
#define ACB_HYPGEOM_KIND_BESSEL_K 6

/* Series and asymptotic expansions, optionally scaled */
int _acb_hypgeom_bessel_asymp_converges(const acb_t nu, const acb_t z, slong prec);
double _acb_hypgeom_bessel_debye_re(double nr, double ni, double wr, double wi);
double _acb_hypgeom_bessel_debye_eta_re(double nr, double ni, double zr, double zi);

/* Re log Gamma(x + iy) to a few digits, for cost and cancellation
   estimates (with regularize set, away from the poles) */
double _acb_hypgeom_d_lgamma_re(double x, double y, int regularize);

/* Estimates in double precision for the choice of algorithm in the
   hypergeometric functions 1F1, U and 2F1 (estimates.c): logarithms of sizes,
   numbers of terms, costs in microseconds. */
typedef struct
{
    double re;
    double im;
}
_acb_hypgeom_cd;

_acb_hypgeom_cd _acb_hypgeom_cd_mul(_acb_hypgeom_cd x, _acb_hypgeom_cd y);
_acb_hypgeom_cd _acb_hypgeom_cd_div(_acb_hypgeom_cd x, _acb_hypgeom_cd y);
_acb_hypgeom_cd _acb_hypgeom_cd_sqrt(_acb_hypgeom_cd x);
_acb_hypgeom_cd _acb_hypgeom_cd_log(_acb_hypgeom_cd x);
double _acb_hypgeom_est_lgamma_re(_acb_hypgeom_cd x);
double _acb_hypgeom_est_lrise(_acb_hypgeom_cd a, double k);
double _acb_hypgeom_est_1f1_peak(double * kpk, _acb_hypgeom_cd a, _acb_hypgeom_cd b, _acb_hypgeom_cd z);
double _acb_hypgeom_est_1f1_terms(_acb_hypgeom_cd a, _acb_hypgeom_cd b, _acb_hypgeom_cd z, double kpk, double target);
void _acb_hypgeom_est_1f1_logabs(double * first, double * second, _acb_hypgeom_cd a, _acb_hypgeom_cd b, _acb_hypgeom_cd z, double cap);
int _acb_hypgeom_est_u_asymp(double * n, double * more, double * cancel, _acb_hypgeom_cd a, _acb_hypgeom_cd b, _acb_hypgeom_cd z, double need, int real);
void _acb_hypgeom_est_u_logabs(double * first, double * second, _acb_hypgeom_cd a, _acb_hypgeom_cd b, _acb_hypgeom_cd z, double cap);
double _acb_hypgeom_est_2f1_peak(double * kpk, _acb_hypgeom_cd a, _acb_hypgeom_cd b, _acb_hypgeom_cd c, double lz);
double _acb_hypgeom_est_2f1_terms(_acb_hypgeom_cd a, _acb_hypgeom_cd b, _acb_hypgeom_cd c, double lz, double kpk, double target);
double _acb_hypgeom_est_2f1_terms_max(_acb_hypgeom_cd a, _acb_hypgeom_cd b, _acb_hypgeom_cd c, double lz, double kpk, double target, double kmax);
void _acb_hypgeom_est_2f1_logabs(double * first, double * second, _acb_hypgeom_cd a, _acb_hypgeom_cd b, _acb_hypgeom_cd c, _acb_hypgeom_cd z, double cap);
/* (the connection formulas of 2F1 at 1/z, 1/(1-z), 1-z and 1-1/z for
   which = 2, ..., 5: the parameters a, b, c of the two series in w,
   log |w|, the logarithms of their factors and of their largest terms
   (at the index kpk), and of the largest term of the formula) */
typedef struct
{
    _acb_hypgeom_cd a[2];
    _acb_hypgeom_cd b[2];
    _acb_hypgeom_cd c[2];
    double lw;
    double pref[2];
    double pk[2];
    double kpk[2];
    double Lam;
    int limit;
}
_acb_hypgeom_est_2f1_conn_struct;

int _acb_hypgeom_est_2f1_conn(_acb_hypgeom_est_2f1_conn_struct * E, int which, _acb_hypgeom_cd a, _acb_hypgeom_cd b, _acb_hypgeom_cd c, _acb_hypgeom_cd z);
double _acb_hypgeom_est_2f1_conn_terms(double * tper, const _acb_hypgeom_est_2f1_conn_struct * E, double level, double kmax);
slong _acb_hypgeom_est_input_accuracy(const acb_t x);
double _acb_hypgeom_est_series_term_cost(double wp, double bits, int real);
double _acb_hypgeom_est_integration_cost(double prec, int real);
int _acb_hypgeom_est_2u_oscillatory(_acb_hypgeom_cd a, _acb_hypgeom_cd b, _acb_hypgeom_cd z);
/* (the same for the two integrals behind 2F1 for Re(z) < 0, with a, b > c) */
int _acb_hypgeom_est_2f1_left_oscillatory(_acb_hypgeom_cd a, _acb_hypgeom_cd b, _acb_hypgeom_cd c,
    _acb_hypgeom_cd z);

/* Enclosures for inexact input from the value v at the midpoints of
   args, for a function that is analytic on the polydisc with radii R_i
   around the midpoints (R_i is ignored for exact arguments) and bounded
   by M there.  If f is given (it evaluates the function at exact
   arguments), the first derivatives are bounded by divided differences
   at a precision of the bits by which M 2^shrink exceeds |v| (at most
   max_prec); shrink is the number of bits by which the R_i are known
   to be smaller than the scales of the function.  Returns 0 if no
   enclosure is obtained (the boxes are too large for the polydisc, M
   is infinite, or v contains zero), 2 if the enclosure is close to the
   range of the function, and 1 if it is valid but M was too crude for
   that; *need_prec is then set to the accuracy of v with which the
   divided differences could be used, if that is the reason. */
int _acb_hypgeom_inexact_enclosure(acb_t res, slong * need_prec, const acb_t v, acb_srcptr args, slong n, const double * R, const mag_t M, slong shrink, void (*f)(acb_t, acb_srcptr, slong, void *), void * ctx, int real, slong max_prec);
void _acb_hypgeom_inexact_polydisc(acb_ptr balls, acb_srcptr args, slong n, const double * R);

/* The driver: f on inexact input from the value at the midpoints and
   a bound majorant(M, balls, stage, ctx) for |f| on the polydisc
   (stage 0: close to |f| and cheap; 1, ..., ncrude: cheap; -1: close to
   |f| and expensive, at exp_ratio times the cost of a value of f);
   cancel is an estimate of the bits that the series for f lose to
   cancellation (for the costs of values at low precision).  With lazy
   set, no more than twice the cost of the value is spent on an
   enclosure that is close to the range of f.  Returns 0, 1 or 2 as
   above. */
#define ACB_HYPGEOM_INEXACT_MAX_ARGS 4
int _acb_hypgeom_inexact(acb_t res, acb_srcptr args, slong n, const double * R,
    void (*f)(acb_t, acb_srcptr, slong, void *),
    void (*majorant)(mag_t, acb_srcptr, int, void *), int ncrude,
    void * ctx, int real, slong prec, slong goal, double cancel, double exp_ratio, int lazy,
    int stage0_cheap);

/* A bound for |pFq(a, b, z)| (p <= q + 1 <= 3, with the factor 1 / k!
   in the terms) on balls from bounds for the absolute values of the
   terms, in double and mag arithmetic; infinite if the series does
   not converge on the balls or needs more than max_terms terms, or
   once the sum exceeds limit (if this is not NULL: a bound of no use
   to the caller, which has a smaller one).  Returns the number of
   terms which were added. */
slong _acb_hypgeom_inexact_terms_bound(mag_t M, acb_srcptr a, slong p, acb_srcptr b, slong q,
    const acb_t z, slong max_terms, const mag_t limit);

/* Whether x is inexact and less accurate than the precision by more
   than 16 bits (more precision does not improve the result then) */
int _acb_hypgeom_inexact_is_coarse(const acb_t x, slong prec);

/* A function f of the inexact argument z (the other arguments are
   exact, in ctx) which is analytic on the box z, from its Taylor
   expansion at the midpoint m:
       f(m + d) = sum_{j < J} c_j d^j + R_J,   |R_J| <= sup |c_J| |d|^J,
   with c_j = f^(j) / j! and the supremum over the box.
   coeff(res, j, x, prec, crude, ctx) sets res to c_j at x: with
   crude = 0, x is the midpoint and the value should be accurate to
   about prec bits; with crude = 1, x is the box (or the midpoint) and
   any enclosure will do, such as the terms summed in ball arithmetic,
   which is too large by their cancellation.  If recur is given,
   recur(res, c_j, c_{j+1}, j, m, prec, ctx) sets res to c_{j+2} at m
   (from the differential equation of f); the coefficients from the
   second on are then computed with it in ball arithmetic, and by coeff
   only where it loses its accuracy.  The order J <= max_order is the
   first for which the remainder with the crude bound is below the
   accuracy goal (or an eighth of the first-order term); *order is set
   to it.  The cost is that of the values c_0 and c_1 at the midpoint
   (at most max_values values with those which replace the
   recurrence), and of J crude values at 64 bits; if the crude values
   show that no order will do, nothing else is computed.  Returns 0 if
   there is no result (no such order, a value which is not finite, or
   f(m) containing zero), and 1 otherwise. */
int _acb_hypgeom_inexact_taylor(acb_t res, slong * order, const acb_t z,
    void (*coeff)(acb_t, slong, const acb_t, slong, int, void *),
    void (*recur)(acb_t, const acb_t, const acb_t, slong, const acb_t, slong, void *),
    void * ctx, int real, slong prec, slong goal, slong max_order, slong max_values);

/* 1F1 by the series or the asymptotic expansion, chosen from the size of z
   only (accurate for parameters of moderate size); U and 2F1 likewise by
   the methods chosen from the location of z (the algorithms of
   acb_hypgeom_u and acb_hypgeom_2f1 without the cost model and the
   integrals; also called by the profile program p-catalog) */
void _acb_hypgeom_m_nointegration(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec);
void _acb_hypgeom_u_nointegration(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec);
void _acb_hypgeom_2f1_nointegration(acb_t res, const acb_t a, const acb_t b, const acb_t c, const acb_t z, int flags, slong prec);
void _acb_hypgeom_bessel_i_nointegration(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec);
void _acb_hypgeom_bessel_k_nointegration(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec);

/* Numerical integration along steepest-descent contours (see sd.h):
   the Bessel function of the given kind (with ACB_HYPGEOM_KIND_BESSEL_JY,
   J and Y in res and res2; otherwise res2 is unused), I and K optionally
   scaled by exp(-+z), and the upper or lower incomplete gamma function.
   The result encloses the function on the input balls.  These return 0
   (leaving res undefined) if the parameters are unsupported or no
   contour is found. */
int _acb_hypgeom_bessel_integration(acb_t res, acb_t res2, int kind, const acb_t nu, const acb_t z, int scaled, slong prec);
int _acb_hypgeom_gamma_incomplete_integration(acb_t res, const acb_t s, const acb_t z, int upper, slong prec);

/* Bounds for sup |F|, sup |dF/dz| and sup |dF/dnu| (F = J, Y, I or K),
   and for sup |F| and sup |dF/ds| (F an incomplete gamma function) over
   the input balls, from the L1 norms of the integrands along the same
   contours; these return 0 on failure (the bounds set to infinity). */
int _acb_hypgeom_bessel_integration_bounds(mag_t M, mag_t Dz, mag_t Dnu, int kind, const acb_t nu, const acb_t z);
int _acb_hypgeom_gamma_incomplete_integration_bounds(mag_t M, mag_t Ds, const acb_t s, const acb_t z, int upper);

/* res (and res2 unless NULL) = the union of f over the parts Im >= 0 and
   Im <= 0 (side = 1, -1) of a ball z crossing the branch cut (-inf, 0) */
int _acb_hypgeom_split_at_cut(acb_t res, acb_t res2, const acb_t nu, const acb_t z, slong prec,
    int (*f)(acb_t, acb_t, const acb_t, const acb_t, int, slong, void *), void * param);

/* The integral of exp(k + lambda t + sum_{j < m} mu_j Log(al_j + be_j t))
   along a polygon with exact vertices (exp_log_integration.c).  Returns
   0 if the integral is not computed. */
int _acb_hypgeom_exp_log_path(acb_t res, const acb_t lambda, acb_srcptr mu, acb_srcptr al,
    acb_srcptr be, const acb_t konst, const acb_t g1, int m, int end0, int end1, acb_srcptr pts,
    slong np, const acb_t c, double logtol, slong prec);

/* The integral of exp(lambda t + sum_j mu_j Log(al_j + be_j t)) over
   [0, 1], or over the ray from 0 in the direction dir (a Gaussian
   integer) if infinite, along a path of steepest descent.  Returns 0 if
   the integral is not computed. */
int _acb_hypgeom_integrate_exp_log(acb_t res, const acb_t dir, const acb_t lambda, acb_srcptr mu,
    acb_srcptr al, acb_srcptr be, int m, int end0, int end1, int infinite, double logscale,
    slong prec);

/* The same over the polygon from 0 by way of the exact number via to 1,
   in place of the segment [0, 1]: the continuation of the integral to
   the other side of a point where a term vanishes, if that point lies
   between the segment and the polygon.  The path is planned with the
   kernel as it is and, if none is found, with the singular points
   moved away from the side of via; in the other order if tilted is
   nonzero (a singular point on the segment or close to it). */
int _acb_hypgeom_integrate_exp_log_via(acb_t res, const acb_t via, const acb_t lambda, acb_srcptr mu,
    acb_srcptr al, acb_srcptr be, int m, int end0, int end1, double logscale, int tilted,
    slong prec);

/* The Euler integrals int_0^1 e^(z t) t^a1 (1 - t)^ba1 dt,
   Gamma(a) U(a, b, z) = int_0^inf e^(-z t) t^a1 (1 + t)^ba1 dt (a1 = a - 1,
   ba1 = b - a - 1; over a ray on which the integral converges) and
   int_0^1 t^b1 (1 - t)^cb1 (1 - z t)^nega dt for complex parameters with
   Re(a1), Re(ba1), Re(b1), Re(cb1) >= 0; these return 0 if the integral
   is not computed. */
int _acb_hypgeom_1f1_integral(acb_t res, const acb_t a1, const acb_t ba1, const acb_t z, slong prec);
int _acb_hypgeom_u_integral(acb_t res, const acb_t a1, const acb_t ba1, const acb_t z, slong prec);
int _acb_hypgeom_2f1_integral(acb_t res, const acb_t b1, const acb_t cb1, const acb_t nega,
    const acb_t z, slong prec);

/* acb_hypgeom_{1f1,u,2f1}_integration without the real-line integral for
   real input (for callers which have tried that already) */
void _acb_hypgeom_1f1_integration_paths(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec);
void _acb_hypgeom_u_integration_paths(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec);
void _acb_hypgeom_2f1_integration_paths(acb_t res, const acb_t a, const acb_t b, const acb_t c, const acb_t z, int regularized, slong prec);

#endif
