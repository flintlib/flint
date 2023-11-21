/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "double_extras.h"
#include "mag.h"

/*
This is a bad implementation the logarithm function,
Any decent libm should compute log() with < 1 ulp
error, but we have no way to know that this is the case
(the IEEE 754 and C language standards make virtually
no guarantees whatsoever about transcendental functions).

So here is a bad implementation which almost certainly
has much worse accuracy than the libm log(), but which
is simple and for which we can give an explicit error
bound. The only precomputed constants are rational
numbers and logarithms of rational numbers, which are
easy to validate.

First write x = 2^n * t with sqrt(2)/2 < t < sqrt(2)
(a few ulps outside this range is fine). Then write
t = m/32 + u with |u| <= 1/64. We have 23 <= m <= 45.
Let v = 32*u/m (note that |v| <= 1/46). Then we have

log(x) = log(2^n*t)
       = n*log(2) + log(t)
       = n*log(2) + log(m/32 + u)
       = n*log(2) + log(m/32 + m/32 * 32*u/m)
       = n*log(2) + log(m/32) + log(1 + 32*u/m)
       = n*log(2) + (log(m/32) + log(1 + v))

We compute log(x) as t1 + (t2 + t3) where t1 ~= n*log(2)
involves two roundings, t2 involves one rounding,
and t3 ~= v - [v^2/2 + ...] involves two roundings for
the multiplication by 1/m, one rounding for subtraction,
plus a much smaller error contribution due to the
error in the approximation of the terms in brackets.

Now consider the outer sum t1 + (t2 + t3). When n != 0,
t1 will always be about twice as large as (t2+t3).
When n = 0, t2 will always be exactly zero or
about twice as large as t3. So there can be no significant
cancellation, and most of the final error comes from
the largest nonzero term.

By inspection, the final error is clearly bounded by
10 ulp (a 3 ulp or maybe even 2 ulp bound could probably
be proved with a careful calculation).

Multiplying the output by 1 +/- 1e-14, we get guaranteed
upper/lower bounds for the logarithm. This factor
is chosen conservatively so that we get correct bounds
even on x87 or if the CPU rounding mode has been changed.

*/

#define LOG_TAB_STEP 32

const double d_log_tab[] = {
    -0.6931471805599453094172, /* log(16/32) */
    -0.6325225587435104668366, /* log(17/32) */
    -0.5753641449035618548784, /* log(18/32) */
    -0.5212969236332860870771, /* log(19/32) */
    -0.4700036292457355536509, /* log(20/32) */
    -0.4212134650763035505856, /* log(21/32) */
    -0.374693449441410693607, /* log(22/32) */
    -0.3302416868705768562794, /* log(23/32) */
    -0.2876820724517809274392, /* log(24/32) */
    -0.2468600779315257978846, /* log(25/32) */
    -0.2076393647782445016154, /* log(26/32) */
    -0.1698990367953974729004, /* log(27/32) */
    -0.1335313926245226231463, /* log(28/32) */
    -0.09844007281325251990289, /* log(29/32) */
    -0.06453852113757117167292, /* log(30/32) */
    -0.031748698314580301157, /* log(31/32) */
    0.0, /* log(32/32) */
    0.03077165866675368837103, /* log(33/32) */
    0.06062462181643484258061, /* log(34/32) */
    0.08961215868968713261995, /* log(35/32) */
    0.1177830356563834545388, /* log(36/32) */
    0.1451820098444978972819, /* log(37/32) */
    0.1718502569266592223401, /* log(38/32) */
    0.1978257433299198803626, /* log(39/32) */
    0.2231435513142097557663, /* log(40/32) */
    0.2478361639045812567806, /* log(41/32) */
    0.2719337154836417588317, /* log(42/32) */
    0.2954642128938358763867, /* log(43/32) */
    0.3184537311185346158102, /* log(44/32) */
    0.3409265869705932103051, /* log(45/32) */
    0.3629054936893684531378, /* log(46/32) */
    0.3844116989103320397348, /* log(47/32) */
};

const double d_log_inverses[] = {
    0.0,
    1.0, /* 1/1 */
    0.5, /* 1/2 */
    0.3333333333333333333333, /* 1/3 */
    0.25, /* 1/4 */
    0.2, /* 1/5 */
    0.1666666666666666666667, /* 1/6 */
    0.1428571428571428571429, /* 1/7 */
    0.125, /* 1/8 */
    0.1111111111111111111111, /* 1/9 */
    0.1, /* 1/10 */
    0.09090909090909090909091, /* 1/11 */
    0.08333333333333333333333, /* 1/12 */
    0.07692307692307692307692, /* 1/13 */
    0.07142857142857142857143, /* 1/14 */
    0.06666666666666666666667, /* 1/15 */
    0.0625, /* 1/16 */
    0.05882352941176470588235, /* 1/17 */
    0.05555555555555555555556, /* 1/18 */
    0.05263157894736842105263, /* 1/19 */
    0.05, /* 1/20 */
    0.04761904761904761904762, /* 1/21 */
    0.04545454545454545454545, /* 1/22 */
    0.0434782608695652173913, /* 1/23 */
    0.04166666666666666666667, /* 1/24 */
    0.04, /* 1/25 */
    0.03846153846153846153846, /* 1/26 */
    0.03703703703703703703704, /* 1/27 */
    0.03571428571428571428571, /* 1/28 */
    0.03448275862068965517241, /* 1/29 */
    0.03333333333333333333333, /* 1/30 */
    0.03225806451612903225806, /* 1/31 */
    0.03125, /* 1/32 */
    0.0303030303030303030303, /* 1/33 */
    0.02941176470588235294118, /* 1/34 */
    0.02857142857142857142857, /* 1/35 */
    0.02777777777777777777778, /* 1/36 */
    0.02702702702702702702703, /* 1/37 */
    0.02631578947368421052632, /* 1/38 */
    0.02564102564102564102564, /* 1/39 */
    0.025, /* 1/40 */
    0.02439024390243902439024, /* 1/41 */
    0.02380952380952380952381, /* 1/42 */
    0.02325581395348837209302, /* 1/43 */
    0.02272727272727272727273, /* 1/44 */
    0.02222222222222222222222, /* 1/45 */
    0.02173913043478260869565, /* 1/46 */
    0.02127659574468085106383, /* 1/47 */
};

double
mag_d_bad_log(double x)
{
    double t, u, v, t1, t2, t3;
    int m, n;

    if (x <= 0.0 || (x-x) != (x-x))
    {
        if (x == 0.0)
            return -D_INF;
        else if (x > 0.0)
            return D_INF;
        else
            return D_NAN;
    }

    if (x < 1 + 1./LOG_TAB_STEP && x > 1 - 1./LOG_TAB_STEP)
    {
        return -d_polyval(d_log_inverses, 12, 1.-x);
    }

    t = frexp(x, &n);

    if (t < 0.7071067811865475244008)
    {
        t *= 2.0;
        n -= 1;
    }

    m = floor(t * LOG_TAB_STEP + 0.5);
    u = t * LOG_TAB_STEP - m;

    /* v = -u / m; */
    v = -u * d_log_inverses[m];

    t1 = n * (-d_log_tab[0]);
    t2 = d_log_tab[m - LOG_TAB_STEP / 2];
    t3 = -d_polyval(d_log_inverses, 11, v);

    return t1 + (t2 + t3);
}

double
mag_d_log_upper_bound(double x)
{
    if (x >= 1.0)
        return mag_d_bad_log(x) * (1 + 1e-14);
    else
        return mag_d_bad_log(x) * (1 - 1e-14);
}

double
mag_d_log_lower_bound(double x)
{
    if (x >= 1.0)
        return mag_d_bad_log(x) * (1 - 1e-14);
    else
        return mag_d_bad_log(x) * (1 + 1e-14);
}

