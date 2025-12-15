/*
    Copyright (C) 2012, 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "ulong_extras.h"
#include "arf.h"
#include "fmpq.h"
#include "acb_poly.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "nfloat.h"

slong
_acb_get_mid_mag(const acb_t z)
{
    slong rm, im;

    rm = arf_abs_bound_lt_2exp_si(arb_midref(acb_realref(z)));
    im = arf_abs_bound_lt_2exp_si(arb_midref(acb_imagref(z)));

    return FLINT_MAX(rm, im);
}

slong
_acb_get_rad_mag(const acb_t z)
{
    slong rm, im;

    /* TODO: write mag function */
    arf_t t;
    arf_init(t);

    arf_set_mag(t, arb_radref(acb_realref(z)));
    rm = arf_abs_bound_lt_2exp_si(t);

    arf_set_mag(t, arb_radref(acb_imagref(z)));
    im = arf_abs_bound_lt_2exp_si(t);

    arf_clear(t);

    return FLINT_MAX(rm, im);
}

/* Compute res[0, ..., n-1] = {i} indexing the upper convex hull of {i, y[i]}. */
static slong convex_hull(slong * res, const double * y, slong len)
{
    slong i, n = 0;

    for (i = 0; i < len; i++)
    {
        while (n >= 2 && (res[n - 2] - res[n - 1]) * (y[i] - y[res[n - 1]])
                <= (i - res[n - 1]) * (y[res[n - 2]] - y[res[n - 1]]))
            n--;

        res[n] = i;
        n++;
    }

    return n;
}

void
_acb_poly_roots_initial_values(acb_ptr roots, acb_srcptr poly, slong deg, slong prec)
{
    double * alog;
    mag_ptr amag;
    mag_t r;
    arf_t ar;
    slong i, j;
    slong * ki, num, m, total;
    double theta;
    acb_t cmid;
    acb_ptr ri;

    amag = _mag_vec_init(deg + 1);
    alog = flint_malloc(sizeof(double) * (deg + 1));
    ki = flint_malloc((deg + 1) * sizeof(slong));
    mag_init(r);
    arf_init(ar);
    acb_init(cmid);  /* shallow only; will not be cleared */

    for (i = 0; i <= deg; i++)
    {
        /* shallow midpoint */
        *arb_midref(acb_realref(cmid)) = *arb_midref(acb_realref(poly + i));
        *arb_midref(acb_imagref(cmid)) = *arb_midref(acb_imagref(poly + i));
        acb_get_mag(amag + i, cmid);

        /* todo: mag_get_d_log2_approx is not very precise; probably
           this does not matter though? */
        alog[i] = mag_get_d_log2_approx(amag + i);
    }

    num = convex_hull(ki, alog, deg + 1);
    total = 0;

    for (i = 1; i < num; i++)
    {
        /* multiplicity */
        m = ki[i] - ki[i - 1];

        /* radius */
        mag_div(r, amag + ki[i - 1], amag + ki[i]);
        mag_root(r, r, m);

        if (mag_is_zero(r))
            mag_set_ui_2exp_si(r, 1, -prec);
        if (mag_is_inf(r))
            mag_set_ui_2exp_si(r, 1, prec);

        arf_set_mag(ar, r);

#if 0
        flint_printf("initial values: radius %{mag} multiplicity %wd\n", r, m);
#endif

        /* Pick m points on the circle with radius r */
        for (j = 0; j < m; j++)
        {
            /* Hack: the angle 2pi * ((j+1)/m + i/deg) distributes points uniformly
               when there are many radii. However, if two radii are
               very close, points can end up overlapping if (j+1)/m + i/deg
               happens to give the same fraction for two different i. Ensure
               that this doesn't happen. */
            double irrational_factor = 0.99071828182845905;

            theta = 6.283185307179586 * (irrational_factor * (j + 1.0) / m + (double) i / deg) + 0.577216;

            ri = roots + total;
            acb_zero(ri);
            arf_set_d(arb_midref(acb_realref(ri)), cos(theta));
            arf_set_d(arb_midref(acb_imagref(ri)), sin(theta));
            arf_mul(arb_midref(acb_realref(ri)), arb_midref(acb_realref(ri)), ar, MAG_BITS, ARF_RND_DOWN);
            arf_mul(arb_midref(acb_imagref(ri)), arb_midref(acb_imagref(ri)), ar, MAG_BITS, ARF_RND_DOWN);
            total++;
        }
    }

    if (total != deg)
        flint_abort();

    _mag_vec_clear(amag, deg + 1);
    flint_free(alog);
    flint_free(ki);
    mag_clear(r);
    arf_clear(ar);
}

#define USE_ABERTH 0

/* todo: this method often gets called with degree 2 or 3 polynomials;
   consider doing something direct there */
int
_acb_poly_find_roots_iter(gr_ptr w, gr_ptr z, gr_srcptr f, gr_srcptr f_prime, slong deg, slong maxiter, gr_ctx_t fp_ctx, gr_ctx_t acb_ctx, slong prec)
{
    slong iter, i;
    slong rootmag, max_rootmag, correction, max_correction;
    slong sz = fp_ctx->sizeof_elem;
    int status = GR_SUCCESS;

    acb_t t;
    acb_init(t);

    for (iter = 0; iter < maxiter; iter++)
    {
        /* todo: support magnitude extraction in the fp_ctx */
        {
            max_rootmag = -WORD_MAX;
            for (i = 0; i < deg; i++)
            {
                status |= gr_set_other(t, GR_ENTRY(z, i, sz), fp_ctx, acb_ctx);
                rootmag = _acb_get_mid_mag(t);
                max_rootmag = FLINT_MAX(rootmag, max_rootmag);
            }
        }

#if USE_ABERTH
            status |= _gr_poly_refine_roots_aberth(w, f, fprime, deg, z, 1, fp_ctx);
#else
            status |= _gr_poly_refine_roots_wdk(w, f, deg, z, 1, fp_ctx);
#endif

        /* read error estimates */
        max_correction = -ARF_PREC_EXACT;
        for (i = 0; i < deg; i++)
        {
            status |= gr_set_other(t, GR_ENTRY(w, i, sz), fp_ctx, acb_ctx);
            correction = _acb_get_mid_mag(t);
            max_correction = FLINT_MAX(correction, max_correction);
        }

        status |= _gr_vec_sub(z, z, w, deg, fp_ctx);

        /* estimate the correction relative to the whole set of roots */
        max_correction -= max_rootmag;
        /* flint_printf("ITER %wd MAX CORRECTION: %wd\n", iter, max_correction); */

        if (max_correction < -prec / 2)
            maxiter = FLINT_MIN(maxiter, iter + 2);
        else if (max_correction < -prec / 3)
            maxiter = FLINT_MIN(maxiter, iter + 3);
        else if (max_correction < -prec / 4)
            maxiter = FLINT_MIN(maxiter, iter + 4);
    }

    acb_clear(t);

    return status;
}

 /* Compute res[0, ..., n-1] = {i} indexing the upper convex hull of {(x[i], y[i])}. */
static slong convex_hull_upper(slong * res, const double * y, const slong * x, slong len)
{
    slong i, n = 0;

    for (i = 0; i < len; i++)
    {
        while (n >= 2 && (x[res[n - 2]] - x[res[n - 1]]) * (y[i] - y[res[n - 1]])
                <= (x[i] - x[res[n - 1]]) * (y[res[n - 2]] - y[res[n - 1]]))
            n--;

        res[n] = i;
        n++;
    }

    return n;
}

/* Computes the line y = ax + b that minimizes max_i |y[i] - (a*x[i] + b)|.
   Returns the maximal vertical distance. 
   Requires x to be sorted. */
static double
_minimax_line(double * res_a, double * res_b, const double * y, const slong * x, slong len)
{
    slong * hU;
    slong nU, pL, pU;
    slong i, j;
    double dyL, dyU;
    slong dxL, dxU;
    double a, b, rL, rU;
    int last_side = 0; /* 0: Lower hull updated, 1: Upper hull updated */

    if (len < 2)
    {
        *res_a = 0.0;
        *res_b = (len == 1) ? y[0] : 0.0;
        return 0.0;
    }

    /* Compute upper hull into hU + 1 so hU[0] is available as a sentinel */
    hU = flint_malloc((len + 1) * sizeof(slong));
    nU = convex_hull_upper(hU + 1, y, x, len);
    /* Set sentinel to 0. Loop condition handles termination. */
    hU[0] = 0;

    /* The lower hull is reduced to the two extreme points and a sentinel at the end */
    slong hL[3] = {0, len-1, 0};

    pL = 0;
    pU = nU;

    /* Initialize edge differences. */
    dyL = y[hL[1]] - y[hL[0]];
    dxL = x[hL[1]] - x[hL[0]];

    /* Upper hull edge (backward) */
    dyU = y[hU[pU]] - y[hU[pU - 1]];
    dxU = x[hU[pU]] - x[hU[pU - 1]];

    /* Loop until the x-coordinate of the lower hull catches up to the upper hull */
    while (x[hL[pL]] < x[hU[pU]])
    {
        /* Compare slopes */
        if (dyL * dxU < dyU * dxL)
        {
            pL++;
            dyL = y[hL[pL + 1]] - y[hL[pL]];
            dxL = x[hL[pL + 1]] - x[hL[pL]];
            last_side = 0;
        }
        else
        {
            pU--;
            dyU = y[hU[pU]] - y[hU[pU - 1]];
            dxU = x[hU[pU]] - x[hU[pU - 1]];
            last_side = 1;
        }
    }

    i = last_side ? hU[pU + 1] : hL[pL - 1];
    j = last_side ? hU[pU]     : hL[pL];

    /* Compute slope a. */
    a = (y[j] - y[i]) / (x[j] - x[i]);

    /* Compute residuals at the support vertices using slope a */
    rL = y[hL[pL]] - a * x[hL[pL]];
    rU = y[hU[pU]] - a * x[hU[pU]];

    /* Compute b and return max distance */
    b = (rU + rL) * 0.5;

    *res_a = a;
    *res_b = b;

    flint_free(hU);

    return (rU - rL) * 0.5;
}

/* Rescale using rotating calipers on the points (i,-log(|poly[i]|)) */
static double
_best_fit_double(arb_t scale, double * cdp, const acb_srcptr poly, slong len)
{
    double *alog;
    double a, b, mdist;
    slong * nonzero;
    slong i, j, prec, prec2;
    mag_t m;
    acb_t cmid, cscale;
    arb_t c, log2;
    
    arb_init(log2);
    arb_init(c);
    acb_init(cmid); /* shallow copies, not cleared */
    acb_init(cscale);
    mag_init(m);
    prec = 53;
    arb_log_ui(log2, 2, prec);
    prec2 = prec + 2*ceil(log(len)/log(2));
    alog = flint_malloc(sizeof(double) * len);
    nonzero = flint_malloc(sizeof(slong) * len);
    j = 0;
    for (i = 0; i < len; i++)
    {
        *arb_midref(acb_realref(cmid)) = *arb_midref(acb_realref(poly + i));
        *arb_midref(acb_imagref(cmid)) = *arb_midref(acb_imagref(poly + i));
        acb_get_mag(m, cmid);
        if( !mag_is_zero(m) ) {
            alog[j] = mag_get_d_log2_approx(m);
            nonzero[j] = i;
            j++;
        }
    }
    
    mdist =  _minimax_line(&a, &b, alog, nonzero, j);  

    if(mdist > 1022) {
        b += mdist - 1022;
    }
    
    arb_set_d(scale, -a);
    arb_mul  (scale, scale, log2, prec);
    arb_exp  (scale, scale, prec);
    arb_set_d(c, -b);
    arb_mul  (c, c, log2, prec);
    arb_exp  (c, c, prec);
    for(i=0; i<len; i++) {
        acb_mul_arb(cscale, poly+i, c, prec);
        cdp[2*i]   = arf_get_d(arb_midref(acb_realref(cscale)), ARF_RND_NEAR);
        cdp[2*i+1] = arf_get_d(arb_midref(acb_imagref(cscale)), ARF_RND_NEAR);
        arb_mul(c, c, scale, prec2);
    }
   
    flint_free(alog);
    flint_free(nonzero);
    mag_clear(m);
    acb_clear(cscale);
    arb_clear(c);
    arb_clear(log2);
    return mdist;
}

double _acb_poly_find_roots_double(acb_ptr roots, acb_srcptr poly, acb_srcptr initial, slong len, slong maxiter, slong prec)
{
    double *cdz, *cdp, *cdi;
    double max_correction;
    arb_t scale;
    mag_t rel_rad, abs_rad;
    slong i;
    
    arb_init(scale);
    mag_init(rel_rad);
    mag_init(abs_rad);
    cdp = flint_malloc(2*len*sizeof(double));
    cdz = flint_malloc(2*(len - 1)*sizeof(double));

    _best_fit_double(scale, cdp, poly, len);
    if(initial != NULL) {
        cdi = cdz;
        for(i=0; i<len-1; i++) {
            acb_div_arb(roots + i, initial + i, scale, prec);
            cdi[2*i]   = arf_get_d(arb_midref(acb_realref(roots + i)), ARF_RND_NEAR);
            cdi[2*i+1] = arf_get_d(arb_midref(acb_imagref(roots + i)), ARF_RND_NEAR);
        }
    } else {
        cdi = NULL;
    }
    max_correction = cd_poly_find_roots(cdz, cdp, cdi, len, maxiter, ldexp(1,-prec));
    slong num_ops = 100*len; //num_ops should be an upper bound on the number of arithmetic operations in a Newton step
    
    /* Gershgorin circle estimate for the radius of the roots */
    mag_set_d(rel_rad, max_correction*(len-2)/(1-FLINT_MIN(0.5, num_ops*0x1p-53)));
    for(i=0; i<len-1; i++) {
        acb_set_d_d(roots + i, cdz[2*i], cdz[2*i+1]);
        acb_get_mag(abs_rad, roots+i);
        mag_mul(abs_rad, abs_rad, rel_rad);
        acb_add_error_mag(roots+i, abs_rad);
        acb_mul_arb(roots + i, roots + i, scale, prec);
    }
    flint_free(cdz);
    flint_free(cdp);
    arb_clear(scale);
    mag_clear(rel_rad);
    mag_clear(abs_rad);
    return max_correction;
}

slong
_acb_poly_find_roots(acb_ptr roots,
    acb_srcptr poly,
    acb_srcptr initial, slong len, slong maxiter, slong prec)
{
    slong i, deg;
    gr_ptr w, z, f, fprime;
    gr_ctx_t acb_ctx, fp_ctx;
    slong sz;
    acb_t t;
    int status = GR_SUCCESS;
    int attempt;

    deg = len - 1;

    if (deg == 0)
    {
        return 0;
    }
    else if (acb_contains_zero(poly + len - 1))
    {
        /* if the leading coefficient contains zero, roots can be anywhere */
        for (i = 0; i < deg; i++)
        {
            arb_zero_pm_inf(acb_realref(roots + i));
            arb_zero_pm_inf(acb_imagref(roots + i));
        }
        return 0;
    }
    else if (deg == 1)
    {
        acb_inv(roots + 0, poly + 1, prec);
        acb_mul(roots + 0, roots + 0, poly + 0, prec);
        acb_neg(roots + 0, roots + 0);
        return 1;
    }

    if (initial != NULL)
        _acb_vec_set(roots, initial, deg);
    else if (prec > 53)
        _acb_poly_roots_initial_values(roots, poly, deg, prec);

    if (maxiter == 0)
        maxiter = 2 * deg + n_sqrt(prec);

    gr_ctx_init_complex_acb(acb_ctx, prec + 64);
    acb_init(t);

    for (attempt = 0; attempt <= 1; attempt++)
    {
        /* First try with nfloat if possible, otherwise fall back to acf */
        if (attempt == 0)
        {
            if (prec <= 53) {
                slong fmaxiter = FLINT_MAX(150,maxiter);
                double max_correction = _acb_poly_find_roots_double(roots, poly, initial, len, fmaxiter, prec);
                status = (max_correction < ldexp(1, -3)) ? GR_SUCCESS : GR_UNABLE ;
                if(status == GR_SUCCESS)
                    break;
            }
            if (nfloat_complex_ctx_init(fp_ctx, prec, 0) != GR_SUCCESS)
                continue;
        }
        else
        {
#if 0
            flint_printf("second try: %wd\n", prec);
#endif
            gr_ctx_init_complex_float_acf(fp_ctx, prec);
        }

        status = GR_SUCCESS;

        sz = fp_ctx->sizeof_elem;
        z = gr_heap_init_vec(4 * deg + 1, fp_ctx);
        w = GR_ENTRY(z, deg, sz);
        fprime = GR_ENTRY(z, 2 * deg, sz);
        f = GR_ENTRY(z, 3 * deg, sz);

        for (i = 0; i <= deg; i++)
            status |= gr_set_other(GR_ENTRY(f, i, sz), poly + i, acb_ctx, fp_ctx);
        for (i = 0; i < deg; i++)
            status |= gr_set_other(GR_ENTRY(z, i, sz), roots + i, acb_ctx, fp_ctx);

#if USE_ABERTH
        status |= _gr_poly_derivative(fprime, f, deg + 1, fp_ctx);
#endif

        if (status == GR_SUCCESS)
            status = _acb_poly_find_roots_iter(w, z, f, fprime, deg, maxiter, fp_ctx, acb_ctx, prec);

        if (status == GR_SUCCESS)
        {
            /* convert back to acb */
            for (i = 0; i < deg; i++)
            {
                GR_MUST_SUCCEED(gr_set_other(roots + i, GR_ENTRY(z, i, sz), fp_ctx, acb_ctx));
                mag_zero(arb_radref(acb_realref(roots + i)));
                mag_zero(arb_radref(acb_imagref(roots + i)));
            }
        }

        gr_heap_clear_vec(z, 4 * deg + 1, fp_ctx);
        gr_ctx_clear(fp_ctx);

        if (status == GR_SUCCESS)
            break;
    }

    acb_clear(t);
    gr_ctx_clear(acb_ctx);

    return _acb_poly_validate_roots(roots, poly, len, prec);
}

slong
acb_poly_find_roots(acb_ptr roots,
    const acb_poly_t poly, acb_srcptr initial,
    slong maxiter, slong prec)
{
    slong len = poly->length;

    if (len == 0)
    {
        flint_throw(FLINT_ERROR, "find_roots: expected a nonzero polynomial");
    }

    return _acb_poly_find_roots(roots, poly->coeffs, initial, len, maxiter, prec);
}
