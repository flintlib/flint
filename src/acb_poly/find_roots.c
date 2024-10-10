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

static slong
_acb_get_mid_mag(const acb_t z)
{
    slong rm, im;

    rm = arf_abs_bound_lt_2exp_si(arb_midref(acb_realref(z)));
    im = arf_abs_bound_lt_2exp_si(arb_midref(acb_imagref(z)));

    return FLINT_MAX(rm, im);
}

#if 0
static slong
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
#endif

/* Compute res[0, ..., n-1] = {i} indexing the convex hull of {i, y[i]}. */
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

static void
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

        /* pick m points on the circle with radius r */
        for (j = 0; j < m; j++)
        {
            theta = 6.283185307179586 * ((j + 1.0) / m + (double) i / deg) + 0.577216;

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
static int
_acb_poly_find_roots_iter(gr_ptr w, gr_ptr z, gr_srcptr f, gr_srcptr FLINT_UNUSED(f_prime), slong deg, slong maxiter, gr_ctx_t fp_ctx, gr_ctx_t acb_ctx, slong prec)
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
                GR_MUST_SUCCEED(gr_set_other(t, GR_ENTRY(z, i, sz), fp_ctx, acb_ctx));
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
            GR_MUST_SUCCEED(gr_set_other(t, GR_ENTRY(w, i, sz), fp_ctx, acb_ctx));
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

    if (initial == NULL)
        _acb_poly_roots_initial_values(roots, poly, deg, prec);
    else
        _acb_vec_set(roots, initial, deg);

    if (maxiter == 0)
        maxiter = 2 * deg + n_sqrt(prec);

    gr_ctx_init_complex_acb(acb_ctx, prec + 64);
    acb_init(t);

    for (attempt = 0; attempt <= 1; attempt++)
    {
        /* First try with nfloat if possible, otherwise fall back to acf */
        if (attempt == 0)
        {
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
