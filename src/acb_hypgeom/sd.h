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
    Steepest-descent contours for integrals of exp(phi(w)) with

        phi(w) = alpha e^w + beta e^(-w) - nu w,

    planned in double precision with rigorous L1 bounds (sd_plan.c), and
    rigorous Gauss-Legendre quadrature along them (sd_integrate.c).  This
    covers the Schlafli integrals of the Bessel functions (alpha = -beta
    = zeta / 2) and the incomplete gamma functions in the logarithmic
    variable (alpha = -1, beta = 0, nu = -s).  Private to acb_hypgeom.
*/

#ifndef ACB_HYPGEOM_SD_H
#define ACB_HYPGEOM_SD_H

#include <math.h>
#include "acb_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SD_MAXSTEPS 3000
#define SD_MAXSADDLES 13       /* including the starting point */
#define SD_MAXHALF (3 * SD_MAXSADDLES)
#define SD_MAXROUTE 10
#define SD_MAXPIECES (2 * SD_MAXROUTE)
#define SD_EPS 2.2e-16
#define SD_XMAX 650.0
#define SD_PI 3.14159265358979323846
#define SD_LN2 0.69314718055994530942

/* complex doubles */
typedef struct { double re, im; } cd;

static inline cd cd_make(double re, double im) { cd z; z.re = re; z.im = im; return z; }
static inline cd cd_add(cd a, cd b) { return cd_make(a.re + b.re, a.im + b.im); }
static inline cd cd_sub(cd a, cd b) { return cd_make(a.re - b.re, a.im - b.im); }
static inline cd cd_mul(cd a, cd b) { return cd_make(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re); }
static inline cd cd_scale(cd a, double s) { return cd_make(a.re * s, a.im * s); }
static inline cd cd_neg(cd a) { return cd_make(-a.re, -a.im); }
static inline cd cd_conj(cd a) { return cd_make(a.re, -a.im); }
static inline double cd_abs(cd a) { return hypot(a.re, a.im); }
static inline double cd_arg(cd a) { return atan2(a.im, a.re); }
static inline cd cd_expi(double t) { return cd_make(cos(t), sin(t)); }

/* The problem: the midpoints of alpha, beta, nu as doubles and radii
   covering the input balls (and the rounding of the midpoints) */
typedef struct
{
    cd al, be, nu;
    double ral, rbe, rnu;
    double aal, abe, anu;    /* absolute values */
    double arg_al, arg_mbe;  /* arg alpha, arg(-beta) */
    cd lin_dir;              /* direction of the linear valley ... */
    double lin_mu;           /* ... and its rate of decay (<= 0: none) */
    int has_start;           /* a finite starting point: */
    cd start;                /* its double approximation, */
    double start_rad;        /* and the radius around it of a disk
                                containing the starting point */
    double drop;             /* trace the descent paths finely until
                                exp(Re phi - top) < e^-drop */
    double tailtol;          /* ... and until the tails are below e^-tailtol
                                times the L1 norm of the polygon */
    double maxdrop;          /* the largest drop of Re phi per step */
}
sd_problem;

/* The endpoints of a contour: the valleys

       SD_PLUS(k):   Re w -> +inf, Im w near pi - arg(alpha) + 2 pi k
       SD_MINUS(k):  Re w -> -inf, Im w near arg(-beta) + 2 pi k   (beta != 0)
       SD_LINEAR:    w -> inf in a direction e^(i theta) along which
                     Re(nu w) -> +inf while alpha e^w + beta e^(-w) stays
                     bounded (cos theta = 0, or cos theta < 0 if beta = 0)

   or the finite starting point of the problem; a descent path can
   also end at another saddle (SD_SADDLE, internally) */
#define SD_PLUS 1
#define SD_MINUS 2
#define SD_LINEAR 3
#define SD_SADDLE 4
#define SD_START 5

typedef struct { int type; int k; } sd_end;

static inline sd_end sd_end_make(int type, int k) { sd_end v; v.type = type; v.k = k; return v; }

/* A descent path from a saddle (or from the starting point), as a
   polygon, with L1 bounds for the weights

       1,   e^x + [beta != 0] e^(-x),   |w|       (x = Re w)

   as logarithms relative to top = Re phi at the saddle */
typedef struct
{
    int saddle;
    sd_end label;            /* where it ends: a valley whose tail bound
                                was certified, or a saddle; type 0 if
                                the path could not be traced */
    double logL[3];          /* L1 bounds with the tail, per weight */
    double logtail[3];       /* the tail bounds alone */
    cd * pts;                /* the vertices, from the saddle */
    int npts;
}
sd_half;

typedef struct
{
    cd w;
    int kind;                /* 2: simple, 3: degenerate (cubic), 1: starting point */
    double width;            /* initial step length */
    cd c;                    /* phi(w) */
    int nhalf;
    int half[3];
    int traced;
}
sd_saddle;

typedef struct
{
    sd_saddle saddles[SD_MAXSADDLES];
    sd_half halves[SD_MAXHALF];
    int nsaddles, nhalf;
}
sd_work;

/* A contour between two ends: pieces to be integrated in the order of
   their vertices and added with signs.  The integrand on a piece is
   scaled by exp(-c) for an exact complex double c (the value of phi at
   the saddle, where |exp(phi - c)| <= 1 or so along the piece); the
   tails beyond the ends of each polygon are bounded by exp(logtail)
   relative to exp(Re c).  Pieces flagged with geom are to be divided
   by 1 - q, q = exp(-2 pi i sigma nu), where sigma is the sign of the
   imaginary part of the direction of the linear valley.  The polygons
   belong to the sd_work or to the contour. */
typedef struct
{
    const cd * pts;
    int npts;
    cd c;
    int sign;
    int geom;                /* to be divided by 1 - q (see _sd_plan) */
    double logtail[2];
    int mid;                 /* for the two descent paths of a simple saddle
                                joined into one piece, from one valley to
                                another: the index of the saddle among the
                                vertices; otherwise -1 */
}
sd_piece;

typedef struct
{
    sd_piece piece[SD_MAXPIECES];
    int npieces;
    double logL[3];          /* log of the L1 bounds (absolute) */
    cd * own[SD_MAXPIECES];  /* polygons belonging to the contour */
    int nown;
}
sd_contour;

void _sd_phi(cd * phi, cd * d1, cd * d2, const sd_problem * P, cd w);
double _sd_margin(const sd_problem * P, cd w);
double _sd_panel(double g[2], const sd_problem * P, cd wa, cd wb, double top);
double _sd_logaddexp(double a, double b);
void _sd_mag_set_exp_d(mag_t res, double L);

/* nu and alpha are given; beta = -alpha if sym, otherwise beta = 0.
   Returns 0 if the midpoints are too large or too small for the
   double precision planning. */
int _sd_problem_init(sd_problem * P, const acb_t alpha, int sym, const acb_t nu,
    const acb_t start, double drop, double tailtol, double maxdrop);

void _sd_work_init(sd_work * W, const sd_problem * P);
void _sd_work_clear(sd_work * W);

/* A contour from the end a to the end b: a steepest-descent route, or
   if there is none and fixed is set, a polygon of straight segments.
   Returns 0 if none is found. */
int _sd_plan(sd_contour * C, sd_work * W, const sd_problem * P, sd_end a, sd_end b, int fixed);
void _sd_contour_clear(sd_contour * C);

/* res = the integral of exp(phi(w)) along the contour, enclosing the
   integral for all parameters in the balls alpha, nu (and starting
   point in the ball start, or NULL); the contour must have been
   planned for the problem P initialized from the same parameters with
   drop and tailtol at least prec log(2) + 12.  Returns 0 on failure. */
int _sd_integrate(acb_t res, const sd_contour * C, const sd_problem * P,
    const acb_t alpha, int sym, const acb_t nu, const acb_t start, slong prec);

/* The same integral by the saddle-point expansions (sd_expand.c), for
   the symmetric kernel (sym = 1) without a starting point, exact alpha
   and nu, and a contour consisting of full descent paths through simple
   saddles which are far enough from the other saddles for the
   expansions to reach the precision.  Returns 0 without computing
   anything if this is not the case, or if the estimated cost exceeds
   maxratio times the estimated cost of the quadrature. */
int _sd_expand(acb_t res, const sd_contour * C, const sd_problem * P,
    const acb_t alpha, int sym, const acb_t nu, slong prec, double maxratio);

/* Estimates for choosing a method (microseconds on a typical machine):
   the cost of the expansion at a saddle with |phi''| = absA and
   |nu / phi''| = q at the given precision (infinite if it does not
   reach the precision; real: with real variables, which cost less),
   and the cost of the quadrature along the two descent paths of a
   saddle with the given number of bits in the arithmetic. */
double _sd_expand_cost(double absA, double q, int real, slong prec);
double _sd_quad_cost(double bits);

#ifdef __cplusplus
}
#endif

#endif
