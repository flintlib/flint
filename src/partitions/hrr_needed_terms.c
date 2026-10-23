/*
    Copyright (C) 2013, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "partitions.h"
#include "partitions/impl.h"

#define INV_LOG2 (1.44269504088896340735992468 + 1e-12)
#define HRR_A (1.1143183348516376904 + 1e-12)  /* 44*pi^2/(225*sqrt(3)) */
#define HRR_B (0.0592384391754448833 + 1e-12)  /* pi*sqrt(2)/75 */
#define HRR_C (2.5650996603237281911 + 1e-12)  /* pi*sqrt(2/3) */

static double
partitions_remainder_bound(double n, double terms)
{
    return HRR_A/sqrt(terms)
            + HRR_B*sqrt(terms/(n-1)) * sinh(HRR_C * sqrt(n)/terms);
}

/* Crude upper bound, sufficient to estimate the precision */
static double
log_sinh(double x)
{
    if (x > 4)
        return x;
    else
        return log(x) + x*x*(1/6.);
}

static double
partitions_remainder_bound_log2(double n, double N)
{
    double t1, t2;

    t1 = log(HRR_A) - 0.5*log(N);
    t2 = log(HRR_B) + 0.5*(log(N) - log(n-1)) + log_sinh(HRR_C * sqrt(n)/N);

    return (FLINT_MAX(t1, t2) + 1) * INV_LOG2;
}

/* The smallest N with the remainder bound below 0.4, found first through
   the crude log2 bound.  Both bounds decrease in N, so the first N where
   the crude one is at most 10 is found by doubling and bisection (the
   result of the linear scan, without its N steps: 1.7e9 at n = 1e20). */
slong
partitions_hrr_needed_terms(double n)
{
    slong N, lo, hi;

    for (hi = 1; partitions_remainder_bound_log2(n, hi) > 10; hi *= 2)
        ;
    lo = hi / 2 + 1;
    if (hi == 1)
        lo = 1;
    while (lo < hi)
    {
        slong mid = lo + (hi - lo) / 2;
        if (partitions_remainder_bound_log2(n, mid) > 10)
            lo = mid + 1;
        else
            hi = mid;
    }
    for (N = lo; partitions_remainder_bound(n, N) > 0.4; N++)
        ;
    return N;
}
