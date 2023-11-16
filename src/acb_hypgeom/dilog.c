/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

#ifdef __GNUC__
# define atan2 __builtin_atan2
# define fabs __builtin_fabs
# define log __builtin_log
#else
# include <math.h>
#endif

void
acb_hypgeom_dilog(acb_t res, const acb_t z, slong prec)
{
    double a, b, best, mz, mz1, t, u;
    int algorithm;
    slong acc, inprec;

    if (!acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (acb_is_zero(z))
    {
        acb_zero(res);
        return;
    }

    acc = acb_rel_accuracy_bits(z);
    acc = FLINT_MAX(acc, 0);
    acc = FLINT_MIN(acc, prec);
    prec = FLINT_MIN(prec, acc + 30);
    inprec = prec;

    /* first take care of exponents that may overflow doubles */
    if (arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), -20) <= 0 &&
        arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), -20) <= 0)
    {
        acb_hypgeom_dilog_zero(res, z, prec);
        return;
    }

    if (arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 20) >= 0 ||
        arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 20) >= 0)
    {
        acb_hypgeom_dilog_transform(res, z, 1, prec);
        return;
    }

    prec = 1.005 * prec + 5;

    a = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_DOWN);
    b = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_DOWN);

    best = mz = a * a + b * b;
    algorithm = 0;

    /* if |z| > 0.25, consider expanding somewhere other than the origin */
    if (best > 0.25 * 0.25)
    {
        if (1.0 / mz < best)  /* use 1/z */
        {
            best = 1.0 / mz;
            algorithm = 1;
        }

        mz1 = (a - 1.0) * (a - 1.0) + b * b;

        if (mz1 < best)   /* use 1-z */
        {
            best = mz1;
            algorithm = 2;
        }

        if (mz1 > 0.001 && mz / mz1 < best)  /* use z/(z-1) */
        {
            best = mz / mz1;
            algorithm = 3;
        }

        if (mz1 > 0.001 && 1.0 / mz1 < best)  /* use 1/(1-z) */
        {
            best = 1.0 / mz1;
            algorithm = 4;
        }
    }

    /* do we still have |z| > 0.25 after transforming? */
    if (best > 0.25 * 0.25)
    {
        /* use series with bernoulli numbers (if not too many!) */
        if (prec < 10000)
        {
            /* t = |log(a+bi)|^2 / (2 pi)^2 */
            t = log(a * a + b * b);
            u = atan2(b, a);
            t = (t * t + u * u) * 0.02533029591;

            if (prec > 1000)
                t *= 4.0;   /* penalty at high precision */
            else
                t *= 1.1;  /* small penalty... also helps avoid this
                              method at negative reals where the log branch
                              cut enters (todo: combine with 1-z formula?) */
            if (t < best)
            {
                algorithm = 8;
                best = t;
            }
        }
    }

    /* fall back on expanding at another special point
       (this should only happen at high precision, where we will use
       the bit-burst algorithm) */
    if (best > 0.75 * 0.75)
    {
        b = fabs(b);   /* reduce to upper half plane */

        /* expanding at z0 = i: effective radius |z-i|/sqrt(2) */
        t = ((b - 1.0) * (b - 1.0) + a * a) * 0.5;
        if (t < best)
        {
            best = t;
            algorithm = 5;
        }

        /* expanding at z0 = (1+i)/2: effective radius |z-(1+i)/2|*sqrt(2) */
        t = 1.0 + 2.0 * (a * (a - 1.0) + b * (b - 1.0));
        if (t < best)
        {
            best = t;
            algorithm = 6;
        }

        /* expanding at z0 = 1+i: effective radius |z-(1+i)| */
        t = 2.0 + (a - 2.0) * a + (b - 2.0) * b;
        if (t < best)
        {
            best = t;
            algorithm = 7;
        }
    }

    if (algorithm == 0)
        acb_hypgeom_dilog_zero(res, z, prec);
    else if (algorithm >= 1 && algorithm <= 7)
        acb_hypgeom_dilog_transform(res, z, algorithm, prec);
    else /* (algorithm == 8) */
        acb_hypgeom_dilog_bernoulli(res, z, prec);

    acb_set_round(res, res, inprec);
}

