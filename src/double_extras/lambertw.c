/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "double_extras.h"

#define POLY(p, x) d_polyval((p), sizeof(p) / sizeof(double), (x))

static const double pol1[4] = {
    0.2278634396856248853716, 0.6685854654191353381433,
    0.4670475452404395343887, 0.061184972065242761167 };

static const double pol2[5] = {
    0.2278636537503804204913, 0.8964421845409468074626,
    1.0217927151592500702475, 0.34513102625055769873401,
    0.020801230123523916719604 };

static const double pol3[6] = {
    0.00005767860320327097931, 0.029896654795890461899563,
    0.0378739044968912982585405, 0.00971957088414193124615358,
    0.000488576886695502361566636, 1.150549466178344373015667e-6 };

static const double pol4[5] = {
    0.030306172539339585635388, 0.066596680780796068408204,
    0.035483738872057375987452, 0.00506436278851840340711316,
    0.0001465263028844943142786722 };

static const double pol5[6] = {
    0.00048233868073637531461, 0.004268700087824343609188,
    0.00127714949974214706149789, 0.0000799706171559085390983949,
    1.186347211803672341928371e-6, 2.943454067276155504308283e-9 };

static const double pol6[6] = {
    0.00553288881087242781512, 0.0043904877060733941697614,
    0.00069354549834088964895342, 0.0000288257440032545960408328,
    3.01054066921000066105342e-7, 4.94316029290773314755549e-10 };

static const double pol7[4] = {
    -0.93011683587619427070, -2.9702322028603227386,
    -2.0759083419960793148, -0.042485660005713612806 };

static const double pol8[4] = {
    0.93011683587619458392, 4.3654074566738568022,
    6.1437079650412473506, 2.4613195056093927345 };

static const double pol9[11] = {
    -1.0000000000000000000, 2.3316439815971242034,
    -1.8121878856393634902, 1.9366311144923597554,
    -2.3535512018816145168, 3.0668589010506319129,
    -4.1753356002581771389, 5.8580237298747741488,
    -8.4010322175239773710, 12.250753501314460424,
    -18.100697012472442755 };

static const double pol10[6] = {
    -5.1972986075163593071, -37.478686466672907613,
    -96.155193004929291698, -102.23856988136744607,
    -37.181958033133170210, -0.48504976999675644134 };

static const double pol11[6] = {
    5.1972986074950082685, 45.274634378414741754, 150.20768172029114131,
    233.88699813222871981, 167.13313463159765859, 42.171248374042409414 };


/* avoid overflows in the formula when x is close to 2^EMAX */
#define RESCALE 1.1102230246251565404e-16

static double
halley(double x, double w)
{
    double t, u, v;

    /* exp() does not overflow, since w is an underestimate
        when the asymptotic series is used */
    t = exp(w) * RESCALE;
    u = 2*w + 2;
    v = w*t - x * RESCALE;
    t = w - u*v / (u*t*(w+1) - (w+2)*v);

    return t;
}

/* this should be exactly 6627126856707895 * 2^(-54) ~=
   0.36787944117144228, which
   is the most negative double in the domain */
#define ONE_OVER_E ldexp(6627126856707895.0, -54)

/* difference from -1/e */
#define CORRECTION 4.3082397558469466e-17

double
d_lambertw(double x)
{
    double t, u, w;

    if (x == 0.0 || x != x || x == D_INF)
        return x;

    if (x < 0.0)
    {
        /* complex result */
        if (x < -ONE_OVER_E)
            return D_NAN;
        /* close to zero */
        else if (x > -1e-9)
            return x - x * x;
        /* close to the singularity at -1/e */
        else if (x + ONE_OVER_E < 0.0003)
            return POLY(pol9, sqrt((x + ONE_OVER_E) + CORRECTION));

        /* otherwise get initial value for Halley iteration */
        if (x + ONE_OVER_E < 0.04)
            w = POLY(pol9, sqrt((x + ONE_OVER_E) + CORRECTION));
        else
            w = x * (1.0 + x * POLY(pol10, x) / POLY(pol11, x));
    }
    else
    {
        /* close to zero */
        if (x <= 0.03125)
        {
            if (x < 1e-9)
                return x - x * x;
            else
                return x * (1.0 + x * POLY(pol7, x) / POLY(pol8, x));
        }

        /* get initial value for Halley iteration  */
        if (x <= 1.0)
            w = x * POLY(pol1, x) / POLY(pol2, x);
        else if (x <= 6.0)
            w = POLY(pol3, x) / POLY(pol4, x);
        else if (x <= 40.0)
            w = POLY(pol5, x) / POLY(pol6, x);
        else
        {
            /* asymptotic series */
            t = log(x);
            u = log(t);
            w = (2*t*t*t - 2*(1+(t-1)*t)*u + u*u)/(2*t*t);
            /* one extra refinement */
            if (x < 1e15)
                w = halley(x, w);
        }
    }

    return halley(x, w);
}
