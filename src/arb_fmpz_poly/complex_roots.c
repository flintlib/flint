/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "profiler.h"
#include "fmpz_poly.h"
#include "arb_poly.h"
#include "acb_poly.h"
#include "arb_fmpz_poly.h"

static int check_accuracy(acb_ptr vec, slong len, slong prec)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (acb_rel_accuracy_bits(vec + i) < prec)
            return 0;
    }

    return 1;
}

static int check_isolation(acb_srcptr roots, slong len)
{
    slong i, j;

    for (i = 0; i < len; i++)
    {
        if (arf_sgn(arb_midref(acb_imagref(roots + i))) >= 0)
        {
            for (j = i + 1; j < len; j++)
            {
                if (arf_sgn(arb_midref(acb_imagref(roots + j))) >= 0)
                {
                    if (acb_overlaps(roots + i, roots + j))
                        return 0;
                }
            }
        }
    }

    return 1;
}


/* Whether the sorted vector contains two entries x, y with
   |x - y| < max(|x|, |y|) / (16 n). */
static int
_arb_vec_has_close_pair(arb_srcptr x, slong len, slong n)
{
    arb_t d, e, m;
    slong i;
    int res = 0;

    arb_init(d);
    arb_init(e);
    arb_init(m);

    for (i = 0; i + 1 < len && !res; i++)
    {
        arb_sub(d, x + i + 1, x + i, 30);
        arb_abs(d, d);
        arb_mul_ui(d, d, 16 * n, 30);
        arb_abs(m, x + i);
        arb_abs(e, x + i + 1);
        arb_max(m, m, e, 30);
        if (!arb_gt(d, m))
            res = 1;
    }

    arb_clear(d);
    arb_clear(e);
    arb_clear(m);

    return res;
}

/* Normwise relative accuracy (in bits) of the coefficients of poly. */
static slong
_arb_poly_normwise_accuracy(const arb_poly_t poly)
{
    mag_t mx, mr, t;
    slong i, acc;

    mag_init(mx);
    mag_init(mr);
    mag_init(t);

    for (i = 0; i < poly->length; i++)
    {
        arb_get_mag(t, poly->coeffs + i);
        mag_max(mx, mx, t);
        mag_max(mr, mr, arb_radref(poly->coeffs + i));
    }

    if (mag_is_zero(mr))
        acc = ARF_PREC_EXACT;
    else if (mag_is_zero(mx) || mag_is_inf(mr))
        acc = -ARF_PREC_EXACT;
    else
        acc = (slong) (mag_get_d_log2_approx(mx) - mag_get_d_log2_approx(mr));

    mag_clear(mx);
    mag_clear(mr);
    mag_clear(t);

    return acc;
}

/* Given enclosures real[0], ..., real[nreal - 1] of the real roots of poly
   (to be refined in place as needed), set Q to an enclosure of the
   exact quotient poly / prod (x - real[i]), whose roots are the nonreal
   roots of poly, with a normwise accuracy of about prec bits.
   The precision needed to compensate for cancellation is tracked in *loss.
   Returns 0 if this appears to be hopeless. */
static int
_arb_fmpz_poly_deflate_real_roots(arb_poly_t Q, const fmpz_poly_t poly,
    arb_ptr real, slong nreal, slong prec, slong * loss, slong max_loss)
{
    arb_poly_t F, D, R;
    slong i, dprec, acc;
    int success = 0;

    arb_poly_init(F);
    arb_poly_init(D);
    arb_poly_init(R);

    for (;;)
    {
        dprec = prec + *loss + 10;

        for (i = 0; i < nreal; i++)
            if (arb_rel_accuracy_bits(real + i) < dprec + 10)
                arb_fmpz_poly_refine_root_arb(real + i, poly, real + i, dprec + 10);

        arb_poly_set_fmpz_poly(F, poly, dprec);
        arb_poly_product_roots(D, real, nreal, dprec);
        arb_poly_divrem(Q, R, F, D, dprec);

        acc = _arb_poly_normwise_accuracy(Q);

        if (acc >= prec)
        {
            success = 1;
            break;
        }

        *loss = FLINT_MAX(*loss + 10, dprec - acc);

        if (*loss > max_loss)
            break;
    }

    arb_poly_clear(F);
    arb_poly_clear(D);
    arb_poly_clear(R);

    return success;
}

void
arb_fmpz_poly_complex_roots(acb_ptr roots, const fmpz_poly_t poly, int flags, slong target_prec)
{
    slong i, j, prec, deg, deg_deflated, isolated, maxiter, deflation;
    slong initial_prec, num_real, nreal_deflated, loss, max_prec;
    acb_poly_t cpoly, cpoly_deflated;
    arb_poly_t Q;
    arb_ptr real_deflated;
    timeit_t timer;
    fmpz_poly_t poly_deflated;
    acb_ptr roots_deflated;
    int removed_zero;

    if (fmpz_poly_degree(poly) < 1)
        return;

    initial_prec = 53;

    fmpz_poly_init(poly_deflated);
    acb_poly_init(cpoly);
    acb_poly_init(cpoly_deflated);

    /* try to write poly as poly_deflated(x^deflation), possibly multiplied by x */
    removed_zero = fmpz_is_zero(poly->coeffs);
    if (removed_zero)
        fmpz_poly_shift_right(poly_deflated, poly, 1);
    else
        fmpz_poly_set(poly_deflated, poly);
    deflation = arb_fmpz_poly_deflation(poly_deflated);
    arb_fmpz_poly_deflate(poly_deflated, poly_deflated, deflation);

    deg = fmpz_poly_degree(poly);
    deg_deflated = fmpz_poly_degree(poly_deflated);

    if (flags & ARB_FMPZ_POLY_ROOTS_VERBOSE)
    {
        flint_printf("searching for %wd roots, %wd deflated\n", deg, deg_deflated);
    }

    /* we only need deg_deflated entries, but the remainder will be useful
       as scratch space */
    roots_deflated = _acb_vec_init(deg);

    arb_poly_init(Q);
    real_deflated = _arb_vec_init(deg_deflated);

    /* Isolate the real roots first. Real root isolation is usually much
       cheaper than complex root isolation. If all roots are real, we are
       done. Otherwise, we divide out the real roots and only need to
       compute the nonreal roots numerically. This also removes
       any clusters of real roots which would slow down the Durand-Kerner
       iteration. */
    nreal_deflated = 0;
    if (deg_deflated >= 2)
    {
        nreal_deflated = arb_fmpz_poly_real_roots(real_deflated,
            poly_deflated, 0, FLINT_MAX(target_prec, initial_prec));

        /* Dividing out a few real roots does not reduce the work much,
           and makes the polynomial dense and inexact; it is only
           worthwhile if the real roots are clustered. */
        if (nreal_deflated != deg_deflated && 4 * nreal_deflated < deg_deflated &&
            !_arb_vec_has_close_pair(real_deflated, nreal_deflated, deg_deflated))
        {
            nreal_deflated = 0;
        }
    }

    loss = 0;
    max_prec = 16 * FLINT_MAX(target_prec, 64) + 4 * deg_deflated;

    for (prec = initial_prec; ; prec *= 2)
    {
        if (prec == 106)
            prec = 128;

        if (nreal_deflated != 0 && nreal_deflated != deg_deflated &&
                prec > max_prec)
        {
            /* Give up using the real roots and start over (should
               not happen). */
            if (flags & ARB_FMPZ_POLY_ROOTS_VERBOSE)
                flint_printf("giving up on real roots\n");
            nreal_deflated = 0;
            prec = initial_prec;
        }

        maxiter = FLINT_MIN(4 * deg_deflated + 64, prec);

        /* don't reuse the roots computed with double in case of failure */
        int new_initial = (prec == initial_prec) || (prec  == 128);

        if (flags & ARB_FMPZ_POLY_ROOTS_VERBOSE)
        {
            flint_printf("prec=%wd: ", prec);
            timeit_start(timer);
        }

        if (nreal_deflated == deg_deflated)
        {
            for (i = 0; i < deg_deflated; i++)
                acb_set_arb(roots_deflated + i, real_deflated + i);
            isolated = deg_deflated;
        }
        else if (nreal_deflated != 0)
        {
            slong m = deg_deflated - nreal_deflated;

            isolated = 0;

            if (_arb_fmpz_poly_deflate_real_roots(Q, poly_deflated,
                real_deflated, nreal_deflated, prec, &loss, max_prec))
            {
                acb_poly_set_arb_poly(cpoly_deflated, Q);
                maxiter = FLINT_MIN(4 * m + 64, prec);

                isolated = acb_poly_find_roots(roots_deflated + nreal_deflated,
                    cpoly_deflated, new_initial ? NULL : roots_deflated + nreal_deflated,
                    maxiter, prec);

                /* the roots of Q are known to be nonreal */
                for (i = 0; i < isolated; i++)
                    if (arb_contains_zero(acb_imagref(roots_deflated + nreal_deflated + i)))
                        isolated = 0;

                if (isolated == m)
                {
                    for (i = 0; i < nreal_deflated; i++)
                        acb_set_arb(roots_deflated + i, real_deflated + i);
                    isolated = deg_deflated;
                }
                else
                    isolated = 0;
            }
        }
        else
        {
            acb_poly_set_fmpz_poly(cpoly_deflated, poly_deflated, prec);
            isolated = acb_poly_find_roots(roots_deflated, cpoly_deflated,
                new_initial ? NULL : roots_deflated, maxiter, prec);
        }

        if (flags & ARB_FMPZ_POLY_ROOTS_VERBOSE)
        {
            timeit_stop(timer);
            /* (the format of this line is checked by dev/check_examples.sh) */
            if (nreal_deflated != 0)
                flint_printf("%wd isolated roots (%wd real) | ", isolated, nreal_deflated);
            else
                flint_printf("%wd isolated roots | ", isolated);
            timeit_print(timer, 1);
        }

        if (isolated == deg_deflated)
        {
            if (!check_accuracy(roots_deflated, deg_deflated, target_prec))
                continue;

            if (deflation == 1)
            {
                _acb_vec_set(roots, roots_deflated, deg_deflated);
            }
            else  /* compute all nth roots */
            {
                acb_t w, w2;

                acb_init(w);
                acb_init(w2);

                acb_unit_root(w, deflation, prec);
                acb_unit_root(w2, 2 * deflation, prec);

                for (i = 0; i < deg_deflated; i++)
                {
                    if (arf_sgn(arb_midref(acb_realref(roots_deflated + i))) > 0)
                    {
                        acb_root_ui(roots + i * deflation,
                                    roots_deflated + i, deflation, prec);
                    }
                    else
                    {
                        acb_neg(roots + i * deflation, roots_deflated + i);
                        acb_root_ui(roots + i * deflation,
                            roots + i * deflation, deflation, prec);
                        acb_mul(roots + i * deflation,
                            roots + i * deflation, w2, prec);
                    }

                    for (j = 1; j < deflation; j++)
                    {
                        acb_mul(roots + i * deflation + j,
                                roots + i * deflation + j - 1, w, prec);
                    }
                }

                acb_clear(w);
                acb_clear(w2);
            }

            /* by assumption that poly is squarefree, must be just one */
            if (removed_zero)
                acb_zero(roots + deg_deflated * deflation);

            if (!check_accuracy(roots, deg, target_prec))
                continue;

            /* The real roots are already known if we did real root
               isolation (and did not deflate). */
            if (nreal_deflated == 0 || deflation != 1)
            {
                acb_poly_set_fmpz_poly(cpoly, poly, prec);

                if (!acb_poly_validate_real_roots(roots, cpoly, prec))
                    continue;
            }

            for (i = 0; i < deg; i++)
            {
                if (arb_contains_zero(acb_imagref(roots + i)))
                    arb_zero(acb_imagref(roots + i));
            }

            if (!check_isolation(roots, deg))
            {
                /* extremely unlikely */
                if (flags & ARB_FMPZ_POLY_ROOTS_VERBOSE)
                    flint_printf("isolation failure!\n");

                continue;
            }

            if (flags & ARB_FMPZ_POLY_ROOTS_VERBOSE)
                flint_printf("done!\n");

            break;
        }
    }

    _acb_vec_sort_pretty(roots, deg);

    /* pair conjugates */
    num_real = 0;
    for (i = 0; i < deg; i++)
        if (acb_is_real(roots + i))
            num_real++;

    if (deg != num_real)
    {
        for (i = num_real, j = 0; i < deg; i++)
        {
            if (arb_is_positive(acb_imagref(roots + i)))
            {
                acb_swap(roots_deflated + j, roots + i);
                j++;
            }
        }

        for (i = 0; i < (deg - num_real) / 2; i++)
        {
            acb_swap(roots + num_real + 2 * i, roots_deflated + i);
            acb_conj(roots + num_real + 2 * i + 1, roots + num_real + 2 * i);
        }
    }

    arb_poly_clear(Q);
    _arb_vec_clear(real_deflated, deg_deflated);
    fmpz_poly_clear(poly_deflated);
    acb_poly_clear(cpoly);
    acb_poly_clear(cpoly_deflated);
    _acb_vec_clear(roots_deflated, deg);
}
