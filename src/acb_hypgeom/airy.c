/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "double_extras.h"
#include "acb_hypgeom.h"

#define LOG2 0.69314718055994530942
#define EXP1 2.7182818284590452354

static const double small_log_tab[] = {
  0.0, 0.0, 0.693147180559945309,
  1.09861228866810969, 1.38629436111989062, 1.60943791243410037,
  1.791759469228055, 1.94591014905531331, 2.07944154167983593,
  2.19722457733621938, 2.30258509299404568, 2.39789527279837054,
  2.48490664978800031, 2.56494935746153674, 2.63905732961525861,
  2.70805020110221007, 2.77258872223978124, 2.83321334405621608,
  2.89037175789616469, 2.94443897916644046, 2.99573227355399099,
  3.044522437723423, 3.09104245335831585, 3.13549421592914969,
  3.17805383034794562, 3.21887582486820075, 3.25809653802148205,
  3.29583686600432907, 3.33220451017520392, 3.36729582998647403,
  3.40119738166215538, 3.43398720448514625, 3.46573590279972655,
  3.49650756146648024, 3.52636052461616139, 3.55534806148941368,
  3.58351893845611, 3.61091791264422444, 3.63758615972638577,
  3.66356164612964643, 3.6888794541139363, 3.7135720667043078,
  3.73766961828336831, 3.76120011569356242, 3.78418963391826116,
  3.80666248977031976, 3.828641396489095, 3.85014760171005859,
  3.87120101090789093, 3.89182029811062661, 3.91202300542814606,
  3.93182563272432577, 3.95124371858142735, 3.97029191355212183,
  3.98898404656427438, 4.00733318523247092, 4.02535169073514923,
  4.04305126783455015, 4.06044301054641934, 4.07753744390571945,
  4.09434456222210068, 4.11087386417331125, 4.12713438504509156,
  4.14313472639153269,
};

static slong
asymp_pick_terms(double prec, double logz)
{
    double logk, log2term, log2term_prev;
    slong k;

    log2term_prev = 0.0;

    for (k = 1; ; k++)
    {
        logk = k < 64 ? small_log_tab[k] : log(k);

        log2term = -1.3257480647361593990 - 0.72134752044448170368*logk +
            k * (-1.8577325401678072259 + 1.4426950408889634074*logk
                 - 2.1640425613334451110*logz);

        if (log2term > log2term_prev - 0.5)
            return -1;

        if (log2term < -prec)
            return k;
    }
}

/*
Accurate estimate of log2(|Ai(z)|) or log2(|Bi(z)|), given z = x + yi.
Should subtract 0.25*log2(|z| + log2(2*sqrt(pi)) from the output.
*/
static double
estimate_airy(double x, double y, int ai)
{
    double r, t, sgn;

    r = x;
    t = x * ((x * x) - 3.0 * (y * y));
    y = y * (3.0 * (x * x) - (y * y));
    x = t * (1.0 / 9.0);
    y = y * (1.0 / 9.0);

    sgn = 1.0;

    if (r > 0.0 && x > 0.0)
    {
        t = sqrt(x * x + y * y) + x;

        if (ai) sgn = -1.0;
    }
    else
    {
        x = -x;

        if (x > 0.0 && x > 1e6 * fabs(y))
            t = y * y / (2.0 * x);
        else
            t = sqrt(x * x + y * y) - x;
    }

    return sgn * sqrt(0.5 * t) * 2.8853900817779268147;
}

/* error propagation based on derivatives */
void
acb_hypgeom_airy_prop(acb_t ai, acb_t aip, acb_t bi, acb_t bip,
    const acb_t z, slong n, int algo, slong prec)
{
    mag_t aib, aipb, bib, bipb, zb, rad;
    acb_t zz;
    int real;

    mag_init(aib);
    mag_init(aipb);
    mag_init(bib);
    mag_init(bipb);
    mag_init(zb);
    mag_init(rad);
    acb_init(zz);

    real = acb_is_real(z);
    arf_set(arb_midref(acb_realref(zz)), arb_midref(acb_realref(z)));
    arf_set(arb_midref(acb_imagref(zz)), arb_midref(acb_imagref(z)));
    mag_hypot(rad, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));
    acb_get_mag(zb, z);

    acb_hypgeom_airy_bound(aib, aipb, bib, bipb, z);
    if (algo == 0)
        acb_hypgeom_airy_direct(ai, aip, bi, bip, zz, n, prec);
    else
        acb_hypgeom_airy_asymp(ai, aip, bi, bip, zz, n, prec);

    if (ai != NULL)
    {
        mag_mul(aipb, aipb, rad);
        if (real)
            arb_add_error_mag(acb_realref(ai), aipb);
        else
            acb_add_error_mag(ai, aipb);
    }

    if (aip != NULL)
    {
        mag_mul(aib, aib, rad);
        mag_mul(aib, aib, zb);  /* |Ai''(z)| = |z Ai(z)| */
        if (real)
            arb_add_error_mag(acb_realref(aip), aib);
        else
            acb_add_error_mag(aip, aib);
    }

    if (bi != NULL)
    {
        mag_mul(bipb, bipb, rad);
        if (real)
            arb_add_error_mag(acb_realref(bi), bipb);
        else
            acb_add_error_mag(bi, bipb);
    }

    if (bip != NULL)
    {
        mag_mul(bib, bib, rad);
        mag_mul(bib, bib, zb);  /* |Bi''(z)| = |z Bi(z)| */
        if (real)
            arb_add_error_mag(acb_realref(bip), bib);
        else
            acb_add_error_mag(bip, bib);
    }

    mag_clear(aib);
    mag_clear(aipb);
    mag_clear(bib);
    mag_clear(bipb);
    mag_clear(zb);
    mag_clear(rad);
    acb_clear(zz);
}

void
acb_hypgeom_airy_direct_prop(acb_t ai, acb_t aip, acb_t bi, acb_t bip,
    const acb_t z, slong n, slong prec)
{
    acb_hypgeom_airy_prop(ai, aip, bi, bip, z, n, 0, prec);
}

void
acb_hypgeom_airy_asymp2(acb_t ai, acb_t aip, acb_t bi, acb_t bip,
    const acb_t z, slong n, slong prec)
{
    /* avoid singularity in asymptotic expansion near 0 */
    if (acb_rel_accuracy_bits(z) > 3)
        acb_hypgeom_airy_asymp(ai, aip, bi, bip, z, n, prec);
    else
        acb_hypgeom_airy_prop(ai, aip, bi, bip, z, n, 1, prec);
}

void
acb_hypgeom_airy(acb_t ai, acb_t aip, acb_t bi, acb_t bip, const acb_t z, slong prec)
{
    arf_srcptr re, im;
    double x, y, t, zmag, z15, term_est, airy_est, abstol;
    slong n, wp;

    if (!acb_is_finite(z))
    {
        if (ai != NULL) acb_indeterminate(ai);
        if (aip != NULL) acb_indeterminate(aip);
        if (bi != NULL) acb_indeterminate(bi);
        if (bip != NULL) acb_indeterminate(bip);
        return;
    }

    re = arb_midref(acb_realref(z));
    im = arb_midref(acb_imagref(z));
    wp = prec * 1.03 + 15;

    /* tiny input -- use direct method and pick n without underflowing */
    if (arf_cmpabs_2exp_si(re, -64) < 0 && arf_cmpabs_2exp_si(im, -64) < 0)
    {
        if (arf_cmpabs_2exp_si(re, -wp) < 0 && arf_cmpabs_2exp_si(im, -wp) < 0)
        {
            n = 1;  /* very tiny input */
        }
        else
        {
            if (arf_cmpabs(re, im) > 0)
                zmag = fmpz_get_d(ARF_EXPREF(re));
            else
                zmag = fmpz_get_d(ARF_EXPREF(im));
            zmag = 3.0 * zmag + 1;
            n = wp / (-zmag) + 1;
        }

        if (acb_is_exact(z))
            acb_hypgeom_airy_direct(ai, aip, bi, bip, z, n, wp);
        else
            acb_hypgeom_airy_direct_prop(ai, aip, bi, bip, z, n, wp);
    }  /* huge input -- use asymptotics and pick n without overflowing */
    else if ((arf_cmpabs_2exp_si(re, 64) > 0 || arf_cmpabs_2exp_si(im, 64) > 0))
    {
        if (arf_cmpabs_2exp_si(re, prec) > 0 || arf_cmpabs_2exp_si(im, prec) > 0)
        {
            n = 1;   /* very huge input */
        }
        else
        {
            x = fmpz_get_d(ARF_EXPREF(re));
            y = fmpz_get_d(ARF_EXPREF(im));
            zmag = (FLINT_MAX(x, y) - 2) * LOG2;
            n = asymp_pick_terms(wp, zmag);
            n = FLINT_MAX(n, 1);
        }

        acb_hypgeom_airy_asymp2(ai, aip, bi, bip, z, n, wp);
    }
    else /* moderate input */
    {
        x = arf_get_d(re, ARF_RND_DOWN);
        y = arf_get_d(im, ARF_RND_DOWN);

        zmag = sqrt(x * x + y * y);
        z15 = zmag * sqrt(zmag);

        if (zmag >= 4.0 && (n = asymp_pick_terms(wp, log(zmag))) != -1)
        {
            acb_hypgeom_airy_asymp2(ai, aip, bi, bip, z, n, wp);
        }
        else if (zmag <= 1.5)
        {
            t = 3 * (wp * LOG2) / (2 * z15 * EXP1);
            t = (wp * LOG2) / (2 * d_lambertw(t));
            n = FLINT_MAX(t + 1, 2);

            if (acb_is_exact(z))
                acb_hypgeom_airy_direct(ai, aip, bi, bip, z, n, wp);
            else
                acb_hypgeom_airy_direct_prop(ai, aip, bi, bip, z, n, wp);
        }
        else
        {
            /* estimate largest term: log2(exp(2(z^3/9)^(1/2))) */
            term_est = 0.96179669392597560491 * z15;

            /* estimate the smaller of Ai and Bi */
            airy_est = estimate_airy(x, y, (ai != NULL || aip != NULL));

            /* estimate absolute tolerance and necessary working precision */
            abstol = airy_est - wp;
            wp = wp + term_est - airy_est;
            wp = FLINT_MAX(wp, 10);

            t = 3 * (-abstol * LOG2) / (2 * z15 * EXP1);
            t = (-abstol * LOG2) / (2 * d_lambertw(t));
            n = FLINT_MAX(t + 1, 2);

            if (acb_is_exact(z))
                acb_hypgeom_airy_direct(ai, aip, bi, bip, z, n, wp);
            else
                acb_hypgeom_airy_direct_prop(ai, aip, bi, bip, z, n, wp);
        }
    }

    if (ai != NULL) acb_set_round(ai, ai, prec);
    if (aip != NULL) acb_set_round(aip, aip, prec);
    if (bi != NULL) acb_set_round(bi, bi, prec);
    if (bip != NULL) acb_set_round(bip, bip, prec);
}

