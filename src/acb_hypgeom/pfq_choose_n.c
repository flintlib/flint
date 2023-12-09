/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_hypgeom.h"

#ifdef __GNUC__
# define floor __builtin_floor
# define log __builtin_log
#else
# include <math.h>
#endif

#define D_ABS(x) ((x) < 0.0 ? (-(x)) : (x))

int
acb_hypgeom_pfq_choose_n_double(slong * nn,
    const double * are, const double * aim, slong p,
    const double * bre, const double * bim, slong q,
    double log2_z,
    slong n_skip, slong n_min, slong n_max, slong prec)
{
    double increase, term, term_max, accuracy, accuracy_best, t, u;
    double required_decrease;
    slong k, n, n_best;
    int success;

    if (p < q)
        required_decrease = 0.01;
    else if (p == q)
        required_decrease = 0.0001;
    else
        required_decrease = 0.01;

    term = term_max = accuracy = accuracy_best = 0.0;
    success = 0;

    for (n = n_best = n_skip; n < n_max; n++)
    {
        t = 1.0;

        for (k = 0; k < FLINT_MAX(p, q); k++)
        {
            if (k < p)
            {
                u = (are[k]+n-1)*(are[k]+n-1) + (aim[k]*aim[k]);
                t *= D_ABS(u);
            }

            if (k < q)
            {
                u = (bre[k]+n-1)*(bre[k]+n-1) + (bim[k]*bim[k]);
                u = D_ABS(u);

                if (u > 1e-100)
                    t /= u;
            }
        }

        increase = 0.5 * log(t) * 1.4426950408889634074 + log2_z;

        term += increase;
        term_max = FLINT_MAX(term_max, term);
        accuracy = term_max - term;

        if (accuracy > accuracy_best && n >= n_min && increase < -required_decrease)
        {
            n_best = n;
            accuracy_best = accuracy;
        }

        if (accuracy_best > prec + 4)
        {
            success = 1;
            break;
        }
    }

    *nn = n_best;
    return success;
}

slong
acb_hypgeom_pfq_choose_n_max(acb_srcptr a, slong p,
                         acb_srcptr b, slong q, const acb_t z,
                         slong prec, slong n_max)
{
    slong n_skip, n_min, n_terminating, nint;
    slong k, n;
    double log2_z;
    double * are;
    double * aim;
    double * bre;
    double * bim;
    mag_t zmag;
    int success;

    if (acb_is_zero(z) || !acb_is_finite(z))
        return 1;

    mag_init(zmag);

    are = flint_malloc(sizeof(double) * 2 * (p + q));
    aim = are + p;
    bre = aim + p;
    bim = bre + q;

    acb_get_mag(zmag, z);
    log2_z = mag_get_d_log2_approx(zmag);

    n_skip = 1;
    n_min = 1;
    n_terminating = WORD_MAX;

    for (k = 0; k < p; k++)
    {
        are[k] = arf_get_d(arb_midref(acb_realref(a + k)), ARF_RND_DOWN);
        aim[k] = arf_get_d(arb_midref(acb_imagref(a + k)), ARF_RND_DOWN);

        /* If the series is terminating, stop at this n. */
        if (acb_is_int(a + k) && are[k] <= 0.0)
        {
            n_terminating = FLINT_MIN(n_terminating, (slong) (-are[k] + 1));
            n_terminating = FLINT_MAX(n_terminating, 1);
        }
        else if (are[k] <= 0.01 && D_ABS(aim[k]) < 0.01)
        {
            /* If very close to an integer, we don't want to deal with it
               using doubles, so fast forward (todo: could work with the
               log2 of the difference instead when this happens). */
            nint = floor(are[k] + 0.5);
            if (D_ABS(nint - are[k]) < 0.01)
                n_skip = FLINT_MAX(n_skip, 2 - nint);
        }
    }

    n_max = FLINT_MIN(n_max, n_terminating);

    for (k = 0; k < q; k++)
    {
        bre[k] = arf_get_d(arb_midref(acb_realref(b + k)), ARF_RND_DOWN);
        bim[k] = arf_get_d(arb_midref(acb_imagref(b + k)), ARF_RND_DOWN);

        if (bre[k] <= 0.25)
        {
            n_min = FLINT_MAX(n_min, 2 - bre[k]);

            /* Also avoid near-integers here (can even allow exact
               integers when computing regularized hypergeometric functions). */
            if (bre[k] <= 0.01 && D_ABS(bim[k]) < 0.01)
            {
                nint = floor(bre[k] + 0.5);
                if (D_ABS(nint - bre[k]) < 0.01)
                    n_skip = FLINT_MAX(n_skip, 2 - nint);
            }
        }
    }

    success = acb_hypgeom_pfq_choose_n_double(&n, are, aim, p, bre, bim, q,
        log2_z, n_skip, n_min, n_max, prec);

    if (!success)
    {
        if (n_terminating <= n_max)
        {
            n = n_terminating;
        }
        else
        {
            n = FLINT_MAX(n_min, n);
            n = FLINT_MIN(n_max, n);
        }
    }

    flint_free(are);
    mag_clear(zmag);

    return n;
}

slong
acb_hypgeom_pfq_choose_n(acb_srcptr a, slong p,
                                acb_srcptr b, slong q,
                                const acb_t z, slong prec)
{
    return acb_hypgeom_pfq_choose_n_max(a, p, b, q, z, prec,
        FLINT_MIN(WORD_MAX / 2, 50 + 10 * prec));
}

slong
acb_hypgeom_pfq_series_choose_n(const acb_poly_struct * a, slong p,
                                const acb_poly_struct * b, slong q,
                                const acb_poly_t z, slong len, slong prec)
{
    slong n_skip, n_min, n_max, n_terminating, nint;
    slong k, n;
    double log2_z;
    double * are;
    double * aim;
    double * bre;
    double * bim;
    acb_t t;
    mag_t zmag;
    int success;

    if (z->length == 0) /* todo: check finite */
        return 1;

    mag_init(zmag);
    acb_init(t);

    are = flint_malloc(sizeof(double) * 2 * (p + q));
    aim = are + p;
    bre = aim + p;
    bim = bre + q;

    acb_get_mag(zmag, z->coeffs);
    log2_z = mag_get_d_log2_approx(zmag);

    n_skip = 1;
    n_min = 1;
    n_max = FLINT_MIN(WORD_MAX / 2, 50 + 10 * prec);
    n_terminating = WORD_MAX;

    /* e.g. for exp(x + O(x^100)), make sure that we actually
       run the recurrence up to order 100.
       todo: compute correctly based on degrees of inputs */
    n_min = FLINT_MAX(n_min, len);
    n_max = FLINT_MAX(n_max, n_min);

    for (k = 0; k < p; k++)
    {
        acb_poly_get_coeff_acb(t, a + k, 0);

        are[k] = arf_get_d(arb_midref(acb_realref(t)), ARF_RND_DOWN);
        aim[k] = arf_get_d(arb_midref(acb_imagref(t)), ARF_RND_DOWN);

        /* If the series is terminating, stop at this n. */
        if (acb_is_int(t) && are[k] <= 0.0 && acb_poly_length(a + k) <= 1)
        {
            n_terminating = FLINT_MIN(n_terminating, (slong) (-are[k] + 1));
            n_terminating = FLINT_MAX(n_terminating, 1);
        }
        else if (are[k] <= 0.01 && D_ABS(aim[k]) < 0.01)
        {
            /* If very close to an integer, we don't want to deal with it
               using doubles, so fast forward (todo: could work with the
               log2 of the difference instead when this happens). */
            nint = floor(are[k] + 0.5);
            if (D_ABS(nint - are[k]) < 0.01)
                n_skip = FLINT_MAX(n_skip, 2 - nint);
        }
    }

    n_max = FLINT_MIN(n_max, n_terminating);

    for (k = 0; k < q; k++)
    {
        acb_poly_get_coeff_acb(t, b + k, 0);

        bre[k] = arf_get_d(arb_midref(acb_realref(t)), ARF_RND_DOWN);
        bim[k] = arf_get_d(arb_midref(acb_imagref(t)), ARF_RND_DOWN);

        if (bre[k] <= 0.25)
        {
            n_min = FLINT_MAX(n_min, 2 - bre[k]);

            /* Also avoid near-integers here (can even allow exact
               integers when computing regularized hypergeometric functions). */
            if (bre[k] <= 0.01 && D_ABS(bim[k]) < 0.01)
            {
                nint = floor(bre[k] + 0.5);
                if (D_ABS(nint - bre[k]) < 0.01)
                    n_skip = FLINT_MAX(n_skip, 2 - nint);
            }
        }
    }

    success = acb_hypgeom_pfq_choose_n_double(&n, are, aim, p, bre, bim, q,
        log2_z, n_skip, n_min, n_max, prec);

    if (!success)
    {
        if (n_terminating <= n_max)
        {
            n = n_terminating;
        }
        else
        {
            n = FLINT_MAX(n_min, n);
            n = FLINT_MIN(n_max, n);
        }
    }

    flint_free(are);
    acb_clear(t);
    mag_clear(zmag);

    return n;
}

