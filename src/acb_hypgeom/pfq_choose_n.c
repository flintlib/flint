/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "acb_poly.h"
#include "acb_hypgeom.h"
#include "acb_hypgeom/impl.h"

#define D_ABS(x) ((x) < 0.0 ? (-(x)) : (x))

static int
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
                u = D_ABS(u);

                /* (A factor zero is reached only with a parameter that
                   is a power series whose constant term is a
                   nonpositive integer: the factor is then a multiple
                   of x, and the terms continue without it.) */
                if (u > 1e-100)
                    t *= u;
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

/* (the number of truncation points that are tried before a pole) */
#define F_POLE_CANDIDATES 4

/* The number of terms as in acb_hypgeom_pfq_choose_n_max.  If the bound
   of acb_hypgeom_pfq_bound_factor has been computed for that number,
   it is stored in C and have_bound is set (it costs as much as tens
   of terms when it covers the range up to a pole, and the caller
   needs it again). */
slong
_acb_hypgeom_pfq_choose_n_bound(mag_t C, int * have_bound, acb_srcptr a, slong p,
                         acb_srcptr b, slong q, const acb_t z,
                         slong prec, slong n_max)
{
    slong n_skip, n_min, n_terminating, nint, n_pole = 0;
    slong k, n, n_bound = -1;
    double log2_z;
    double * are;
    double * aim;
    double * bre;
    double * bim;
    mag_t zmag;
    int success;

    *have_bound = 0;

    if (acb_is_zero(z) || !acb_is_finite(z))
        return 1;

    mag_init(zmag);

    are = flint_malloc(sizeof(double) * 2 * (p + q));
    aim = are + p;
    bre = aim + p;
    bim = bre + q;

    acb_get_mag(zmag, z);

    /* (accurately: an error of one unit in log2 |z| would double the
       number of terms for large |z|) */
    if (mag_cmp_2exp_si(zmag, -1000) > 0 && mag_cmp_2exp_si(zmag, 1000) < 0)
        log2_z = log(mag_get_d(zmag)) * 1.4426950408889634074;
    else
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

        /* A pole at a nonpositive integer beyond the last term of a
           terminating series: the terms vanish before the pole, but
           the tail bound is infinite.  Sum all the terms. */
        if (n_terminating == n_max && bre[k] <= 0.25 && 2 - bre[k] >= n_max
            && acb_contains_int(b + k))
            n_min = n_max;

        /* Sum past the negative-integer poles of the Pochhammer
           symbols (unless this is too far out) where the terms grow
           again near the pole to a size that matters; the tail bound
           accounts for the growth otherwise. */
        if (bre[k] <= 0.25 && 2 - bre[k] < n_max)
        {
            n_pole = FLINT_MAX(n_pole, 2 - bre[k]);

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

    if (n_pole > n)
    {
        /* the sizes (log2) of the terms beyond n, up to the pole, and
           at the candidates n + 16, n + 48, n + 112, ... for the
           truncation */
        double term = 0.0, term_max = 0.0, t, u, cterm[F_POLE_CANDIDATES];
        slong m, j, cand[F_POLE_CANDIDATES], nc = 0, next = n;

        for (m = n_skip; m < n_pole; m++)
        {
            t = 1.0;
            for (j = 0; j < FLINT_MAX(p, q); j++)
            {
                if (j < p)
                {
                    u = (are[j]+m-1)*(are[j]+m-1) + (aim[j]*aim[j]);
                    t *= D_ABS(u);
                }
                if (j < q)
                {
                    u = (bre[j]+m-1)*(bre[j]+m-1) + (bim[j]*bim[j]);
                    u = D_ABS(u);
                    if (u > 1e-100)
                        t /= u;
                }
            }
            term += 0.5 * log(t) * 1.4426950408889634074 + log2_z;
            if (m <= n)
                term_max = FLINT_MAX(term_max, term);
            else if (term > term_max - prec - 4)
                break;

            if (m == next && nc < F_POLE_CANDIDATES)
            {
                cand[nc] = m;
                cterm[nc] = term;
                nc++;
                next = n + 16 * ((WORD(1) << nc) - 1);
            }
        }

        /* The terms up to the pole do not matter: the series can be
           truncated before it if the tail bound at the truncation
           point, which covers the range up to the pole with suprema
           of the ratios, is as small as they are (it can be too large
           by many orders of magnitude, or infinite).  A few points
           are tried; the bound decreases with the terms. */
        if (m >= n_pole)
        {
            double lc;

            m = 0;

            for (j = 0; j < nc && m == 0; j++)
            {
                acb_hypgeom_pfq_bound_factor(C, a, p, b, q, z, cand[j]);
                lc = mag_is_finite(C) ? mag_get_d_log2_approx(C) : 1e300;

                if (cterm[j] + lc <= term_max - prec + 4.0)
                {
                    m = n_pole;
                    n_bound = cand[j];

                    if (cand[j] > n)
                    {
                        n_min = FLINT_MAX(n_min, cand[j]);
                        success = acb_hypgeom_pfq_choose_n_double(&n, are, aim, p, bre, bim, q,
                            log2_z, n_skip, n_min, n_max, prec);
                    }
                }
            }
        }

        if (m < n_pole)
        {
            n_min = n_pole;
            success = acb_hypgeom_pfq_choose_n_double(&n, are, aim, p, bre, bim, q,
                log2_z, n_skip, n_min, n_max, prec);
        }
    }

    /* A numerator parameter with negative real part: the terms can
       decrease towards k = -Re(a_i) and grow again beyond, to a second
       maximum.  If the terms beyond n can grow at all, look at their
       sizes, and sum past the last one that matters. */
    if (success && p <= q && n < n_max && n_terminating == WORD_MAX)
    {
        double bound, far = 0.0;
        int scan = 0;

        bound = log2_z;

        for (k = 0; k < q; k++)
        {
            double x0 = bre[k] + n;

            if (x0 <= 0.0)
            {
                bound = -1.0;   /* (the poles are handled above) */
                scan = 0;
                break;
            }

            far += fabs(bre[k]) + fabs(bim[k]);

            if (k < p)
            {
                /* sup |a + j| / |b + j| over j >= n */
                double d = are[k] - bre[k], u, v;
                u = 1.0 + d / x0;
                v = aim[k] / x0;
                u = FLINT_MAX(u * u, 1.0) + v * v;
                bound += 0.5 * log(u) * 1.4426950408889634074;

                if (are[k] < -(double) n)
                {
                    scan = 1;
                    far += 2.0 * (fabs(are[k]) + fabs(aim[k])) + 1.0;
                }
            }
            else
            {
                bound -= 0.5 * log(x0 * x0 + bim[k] * bim[k]) * 1.4426950408889634074;
            }
        }

        if (scan && bound >= 0.0)
        {
            double term = 0.0, term_max = 0.0, t, u, increase;
            slong m, j, last = 0;

            for (m = n_skip; m < n_max; m++)
            {
                t = 1.0;
                for (j = 0; j < FLINT_MAX(p, q); j++)
                {
                    if (j < p)
                    {
                        u = (are[j]+m-1)*(are[j]+m-1) + (aim[j]*aim[j]);
                        t *= D_ABS(u);
                    }
                    if (j < q)
                    {
                        u = (bre[j]+m-1)*(bre[j]+m-1) + (bim[j]*bim[j]);
                        u = D_ABS(u);
                        if (u > 1e-100)
                            t /= u;
                    }
                }

                increase = 0.5 * log(t) * 1.4426950408889634074 + log2_z;
                term += increase;

                if (m <= n)
                    term_max = FLINT_MAX(term_max, term);
                else if (term > term_max - prec - 4)
                    last = m;
                else if (m > far + 2.0 && increase < 0.0)
                    break;      /* decreasing from here on */
            }

            if (last != 0)
            {
                n_min = FLINT_MIN(last + 1, n_max);
                success = acb_hypgeom_pfq_choose_n_double(&n, are, aim, p, bre, bim, q,
                    log2_z, n_skip, n_min, n_max, prec);
            }
        }
    }

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

    *have_bound = (n == n_bound);

    flint_free(are);
    mag_clear(zmag);

    return n;
}

slong
acb_hypgeom_pfq_choose_n_max(acb_srcptr a, slong p,
                         acb_srcptr b, slong q, const acb_t z,
                         slong prec, slong n_max)
{
    mag_t C;
    slong n;
    int have_bound;

    mag_init(C);
    n = _acb_hypgeom_pfq_choose_n_bound(C, &have_bound, a, p, b, q, z, prec, n_max);
    mag_clear(C);

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
        else if (acb_is_int(t) && are[k] <= 0.0 && -are[k] < 1e15)
        {
            /* A nonpositive integer -m with higher coefficients (the
               limit of a parameter, as in the connection formulas of
               2F1): the terms beyond m have a multiple of x in place
               of the factor zero, and are of the size of the other
               factors.  (Fast forwarding to m would give the largest
               number of terms for m beyond that: 4 ms instead of 0.2
               at 333 bits for m = 1000.) */
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
