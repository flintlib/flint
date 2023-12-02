/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

#ifdef __GNUC__
# define ldexp __builtin_ldexp
# define floor __builtin_floor
#else
# include <math.h>
#endif

#define ONE_OVER_PI 0.31830988618379067154
#define PI 3.1415926535897932385

/* We use doubles in a way that keeps the final error <= EPSILON. */
#define EPSILON ldexp(1.0, -30)

/* For |x| <= 2^MAX_EXP, doubles can be used directly for quadrant
   reduction. Note: 2^MAX_EXP must also fit in an int. */
#define MAX_EXP 20

/* Lookup tables, steps of 1/16. */
static const double sin_tab[] = {
  0.0,0.062459317842380198585,0.12467473338522768996,
  0.18640329676226988455,0.24740395925452292960,0.30743851458038085067,
  0.36627252908604756137,0.42367625720393801036,0.47942553860420300027,
  0.53330267353602017333,0.58509727294046215481,0.63460708001526929685,
  0.68163876002333416673,0.72600865526071254966,0.76754350223602703963,
  0.80608110826069299518,0.84147098480789650665,0.87357493516707112023,
  0.90226759409909516292,0.92743691738486767172,0.94898461935558621435,
  0.96682655669618022997,0.98089305702315569609,0.99112919095376166988,
  0.99749498660405443094,0.99996558567824887530,
};

static const double cos_tab[] = {
  1.0,0.99804751070009914963,0.99219766722932905315,
  0.98247331310125525749,0.96891242171064478414,0.95156794804817220215,
  0.93050762191231429115,0.90581368342593642074,0.87758256189037271612,
  0.84592449923106795446,0.81096311950521790219,0.77283494615247154481,
  0.73168886887382088631,0.68768556222050484451,0.64099685816332513036,
  0.59180507509247750546,0.54030230586813971740,0.48668966770196333087,
  0.43117651679866617655,0.37397963082453319229,0.31532236239526866545,
  0.25543376688881169791,0.19454770798898718445,0.13290194445282520566,
  0.070737201667702910088,0.0082962316238583774779,
};

/*  Computes sin(a) and cos(a) and also sets q to the (approximate)
    quadrant of a. The absolute error of the straightforward algorithm
    below can be bounded as the error in the mod pi/2 reduction plus
    the negligible contribution (O(1) * 2^-53) of all other steps.
    Since |a| < 2^20 by assumption, we have 3 bits of margin for the
    final error bound of 2^-30 for the entire algorithm. */
static void
sin_cos(double * sin_a, double * cos_a, int * q, double a)
{
    double as, ac, t, b, bs, bc, v;
    int i, qa;

    *q = qa = floor(a * (2.0 * ONE_OVER_PI));
    a = a - qa * (0.5 * PI);

    if (a < 0.0)
        a = 0.0;
    if (a > 0.5 * PI)
        a = 0.5 * PI;

    i = a * 16.0;

    if (i < 0 || i > 25)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    as = sin_tab[i];
    ac = cos_tab[i];

    b = a - i * (1 / 16.0);
    v = b * b;

    /* Just use the Taylor series. */
    bs = 2.7557319223985890653e-6 * v;
    bs = (-0.0001984126984126984127 + bs) * v;
    bs = (0.0083333333333333333333 + bs) * v;
    bs = (-0.16666666666666666667 + bs) * v;
    bs = (1.0 + bs) * b;

    bc = -2.7557319223985890653e-7 * v;
    bc = (0.000024801587301587301587 + bc) * v;
    bc = (-0.0013888888888888888889 + bc) * v;
    bc = (0.041666666666666666667 + bc) * v;
    bc = (-0.5 + bc) * v;
    bc = 1.0 + bc;

    t = as * bc + ac * bs;
    ac = ac * bc - as * bs;
    as = t;

    if ((qa & 3) == 0)
    {
        *sin_a = as;
        *cos_a = ac;
    }
    else if ((qa & 3) == 1)
    {
        *sin_a = ac;
        *cos_a = -as;
    }
    else if ((qa & 3) == 2)
    {
        *sin_a = -as;
        *cos_a = -ac;
    }
    else
    {
        *sin_a = -ac;
        *cos_a = as;
    }
}

/* FIXME: this is the bottleneck -- an mpn version would be better */
static void
_arb_mod_2pi(arb_t x, slong mag)
{
    arb_t t;
    arf_t q;

    arf_init(q);
    arb_init(t);

    arb_const_pi(t, mag + 53);
    arb_mul_2exp_si(t, t, 1);
    arf_div(q, arb_midref(x), arb_midref(t), mag + 10, ARF_RND_NEAR);
    arf_floor(q, q);
    arb_submul_arf(x, t, q, 53);

    arf_clear(q);
    arb_clear(t);
}

void
_arb_sin_cos_wide(arb_t sinx, arb_t cosx, const arf_t xmid, const mag_t xrad, slong prec)
{
    double m, a, b, r, cos_min, cos_max, sin_min, sin_max;
    double as, ac, bs, bc;
    int i, qa, qb;
    slong mag;

    mag = arf_abs_bound_lt_2exp_si(xmid);

    if (mag > FLINT_MAX(65536, 4 * prec) || mag_cmp_2exp_si(xrad, 3) >= 0)
    {
        if (sinx != NULL)
            arb_zero_pm_one(sinx);
        if (cosx != NULL)
            arb_zero_pm_one(cosx);
        return;
    }
    else if (mag <= MAX_EXP)
    {
        m = arf_get_d(xmid, ARF_RND_DOWN);
        r = mag_get_d(xrad);
    }
    else
    {
        arb_t t;
        arb_init(t);
        arf_set(arb_midref(t), xmid);
        mag_set(arb_radref(t), xrad);
        _arb_mod_2pi(t, mag);

        /* this should not happen */
        if (arf_cmpabs_2exp_si(arb_midref(t), 5) > 0 ||
            mag_cmp_2exp_si(arb_radref(t), 5) > 0)
        {
            flint_printf("unexpected precision loss in sin_cos_wide\n");

            if (sinx != NULL)
                arb_zero_pm_one(sinx);
            if (cosx != NULL)
                arb_zero_pm_one(cosx);
            arb_clear(t);
            return;
        }

        m = arf_get_d(arb_midref(t), ARF_RND_DOWN);
        r = mag_get_d(arb_radref(t));

        arb_clear(t);
    }

    a = m - r;
    b = m + r;

    sin_cos(&as, &ac, &qa, a);
    sin_cos(&bs, &bc, &qb, b);

    sin_min = FLINT_MIN(as, bs);
    sin_max = FLINT_MAX(as, bs);
    cos_min = FLINT_MIN(ac, bc);
    cos_max = FLINT_MAX(ac, bc);

    /* Handle the quadrant crossings. */
    for (i = qa; i < qb; i++)
    {
        if ((i & 3) == 1) cos_min = -1.0;
        if ((i & 3) == 3) cos_max = 1.0;
        if ((i & 3) == 2) sin_min = -1.0;
        if ((i & 3) == 0) sin_max = 1.0;
    }

    if (sinx != NULL)
    {
        a = (sin_max + sin_min) * 0.5;
        r = (sin_max - sin_min) * 0.5 + EPSILON;

        arf_set_d(arb_midref(sinx), a);
        mag_set_d(arb_radref(sinx), r);
        arb_set_round(sinx, sinx, prec);
    }

    if (cosx != NULL)
    {
        a = (cos_max + cos_min) * 0.5;
        r = (cos_max - cos_min) * 0.5 + EPSILON;

        arf_set_d(arb_midref(cosx), a);
        mag_set_d(arb_radref(cosx), r);
        arb_set_round(cosx, cosx, prec);
    }
}

void
arb_sin_cos_wide(arb_t sinx, arb_t cosx, const arb_t x, slong prec)
{
    _arb_sin_cos_wide(sinx, cosx, arb_midref(x), arb_radref(x), prec);
}
