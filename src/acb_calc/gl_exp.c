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
    Building blocks for the rigorous quadrature of integrands of the
    form exp(f(w)) along polygons (steepest-descent contours), shared
    by the integration code of acb_hypgeom and arb_hypgeom:

    - the Gauss-Legendre rules of the cache of acb_calc_integrate as
      arrays of dNb balls (dfloat arithmetic with N <= 4 doubles);

    - the bound for the error of a rule on a chord, given a bound for
      the integrand on a disk containing a Bernstein ellipse;

    - the value of a rule for exp(p(w)) with a polynomial p, in dNb
      arithmetic: where f is given by its Taylor polynomial at a point
      near the chord, with a bound for the rest of the series added to
      the constant coefficient, the exponent is evaluated without the
      cancellation between the large terms of f.
*/

#include "thread_support.h"
#include "acb.h"
#include "acb_calc.h"
#include "acb_calc/impl.h"
#include "dfloat.h"

#define N 1
#define DB(name) d1b_##name
#define DBV(name) _d1b_##name
#define GLE(name) gle_##name##_d1
#include "gl_exp.inc"
#undef N
#undef DB
#undef DBV
#undef GLE

#define N 2
#define DB(name) d2b_##name
#define DBV(name) _d2b_##name
#define GLE(name) gle_##name##_d2
#include "gl_exp.inc"
#undef N
#undef DB
#undef DBV
#undef GLE

#define N 3
#define DB(name) d3b_##name
#define DBV(name) _d3b_##name
#define GLE(name) gle_##name##_d3
#include "gl_exp.inc"
#undef N
#undef DB
#undef DBV
#undef GLE

#define N 4
#define DB(name) d4b_##name
#define DBV(name) _d4b_##name
#define GLE(name) gle_##name##_d4
#include "gl_exp.inc"
#undef N
#undef DB
#undef DBV
#undef GLE

/* The nodes and weights as arb balls come from the cache of
   acb_calc_integrate (the (n + 1) / 2 nonnegative nodes in decreasing
   order, and their weights, possibly with more precision than asked
   for); the cache here holds all the nodes and weights of the rules
   used with dfloat, converted to arrays of dNb balls. */
typedef struct
{
    void * dx[4][ACB_CALC_GL_DFLOAT_STEPS];
    void * dw[4][ACB_CALC_GL_DFLOAT_STEPS];
}
gle_cache_struct;

FLINT_TLS_PREFIX gle_cache_struct * gle_cache = NULL;

static void
gle_cleanup(void)
{
    int N, i;

    if (gle_cache == NULL)
        return;

    for (N = 0; N < 4; N++)
        for (i = 0; i < ACB_CALC_GL_DFLOAT_STEPS; i++)
        {
            flint_free(gle_cache->dx[N][i]);
            flint_free(gle_cache->dw[N][i]);
        }
    flint_free(gle_cache);
    gle_cache = NULL;
}

void
_acb_calc_gl_nodes_dfloat(const void ** x, const void ** w, int N, slong i)
{
    arb_srcptr xa, wa;

    if (gle_cache == NULL)
    {
        gle_cache = flint_calloc(1, sizeof(gle_cache_struct));
        flint_register_cleanup_function(gle_cleanup);
    }

    if (gle_cache->dx[N - 1][i] == NULL)
    {
        slong n = _acb_calc_gl_steps[i];
        void * xx, * ww;

        _acb_calc_gl_nodes(&xa, &wa, i, 53 * N + 30);

        xx = flint_malloc(n * (N + 1) * sizeof(double));
        ww = flint_malloc(n * (N + 1) * sizeof(double));
        if (N == 1) gle_nodes_d1(xx, ww, xa, wa, n);
        else if (N == 2) gle_nodes_d2(xx, ww, xa, wa, n);
        else if (N == 3) gle_nodes_d3(xx, ww, xa, wa, n);
        else gle_nodes_d4(xx, ww, xa, wa, n);
        gle_cache->dx[N - 1][i] = xx;
        gle_cache->dw[N - 1][i] = ww;
    }

    *x = gle_cache->dx[N - 1][i];
    *w = gle_cache->dw[N - 1][i];
}

/* (64 / 15) M ah / ((rho - 1) rho^(2 n - 1)), rho = X + sqrt(X^2 - 1) */
void
_acb_calc_gl_error_bound(mag_t err, const mag_t M, const mag_t ah, double X, slong n)
{
    arb_t rho, t;

    if (!(X > 1.0) || !(X < 1e300))
    {
        mag_inf(err);
        return;
    }

    arb_init(rho);
    arb_init(t);

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

    arb_clear(rho);
    arb_clear(t);
}

void
_acb_calc_gl_exp_poly_dfloat(acb_t res, acb_srcptr coef, slong nc, double are, double aim,
    double bre, double bim, slong i, int N, void * scratch)
{
    const void * x, * w;
    slong n = _acb_calc_gl_steps[i];

    _acb_calc_gl_nodes_dfloat(&x, &w, N, i);
    if (N == 1) gle_exp_poly_d1(res, coef, nc, are, aim, bre, bim, x, w, n, scratch);
    else if (N == 2) gle_exp_poly_d2(res, coef, nc, are, aim, bre, bim, x, w, n, scratch);
    else if (N == 3) gle_exp_poly_d3(res, coef, nc, are, aim, bre, bim, x, w, n, scratch);
    else gle_exp_poly_d4(res, coef, nc, are, aim, bre, bim, x, w, n, scratch);
}
