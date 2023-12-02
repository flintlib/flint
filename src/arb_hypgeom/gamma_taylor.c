/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "double_extras.h"
#include "arb_hypgeom.h"

#define DEBUG 0


const double arb_hypgeom_rgamma_d_tab[128] = {
    1.0,
    0.57721566490153286061,
    -0.65587807152025388108,
    -0.042002635034095235529,
    0.1665386113822914895,
    -0.042197734555544336748,
    -0.0096219715278769735621,
    0.0072189432466630995424,
    -0.0011651675918590651121,
    -0.00021524167411495097282,
    0.00012805028238811618615,
    -0.000020134854780788238656,
    -1.2504934821426706573e-6,
    1.1330272319816958824e-6,
    -2.0563384169776071035e-7,
    6.1160951044814158179e-9,
    5.0020076444692229301e-9,
    -1.1812745704870201446e-9,
    1.0434267116911005105e-10,
    7.782263439905071254e-12,
    -3.6968056186422057082e-12,
    5.100370287454475979e-13,
    -2.0583260535665067832e-14,
    -5.3481225394230179824e-15,
    1.2267786282382607902e-15,
    -1.1812593016974587695e-16,
    1.1866922547516003326e-18,
    1.4123806553180317816e-18,
    -2.2987456844353702066e-19,
    1.7144063219273374334e-20,
    1.3373517304936931149e-22,
    -2.0542335517666727893e-22,
    2.7360300486079998448e-23,
    -1.7323564459105166391e-24,
    -2.3606190244992872873e-26,
    1.8649829417172944307e-26,
    -2.2180956242071972044e-27,
    1.2977819749479936688e-28,
    1.1806974749665284062e-30,
    -1.1245843492770880903e-30,
    1.277085175140866204e-31,
    -7.3914511696151408235e-33,
    1.134750257554215761e-35,
    4.6391346410587220299e-35,
    -5.3473368184391988751e-36,
    3.2079959236133526229e-37,
    -4.4458297365507568821e-39,
    -1.3111745188819887129e-39,
    1.6470333525438138868e-40,
    -1.0562331785035812186e-41,
    2.6784429826430494784e-43,
    2.4247154948517826897e-44,
    -3.736587834535612554e-45,
    2.6283329809401954491e-46,
    -9.2981759953768862996e-48,
    -2.3279424186994705986e-49,
    6.1696208352443874204e-50,
    -4.9282955867709899305e-51,
    2.1835131834145106973e-52,
    -1.2187221891475165553e-54,
    -7.1171088416628746319e-55,
    6.9205040543286892535e-56,
    -3.6764384683566763277e-57,
    8.563098056275654328e-59,
    4.9630454283668443848e-60,
    -7.1542945770816152182e-61,
    4.5517276890885041177e-62,
    -1.6183993053202944344e-63,
    -3.8180434243999502464e-66,
    5.1850524119058482295e-66,
    -4.1671368092239208861e-67,
    1.9162906929373887193e-68,
    -3.8089281324683658733e-70,
    -2.2063861055924121016e-71,
    2.7722310960098954165e-72,
    -1.5987660478100181057e-73,
    5.3197307804174034028e-75,
    -8.0517461416842390432e-78,
    -1.2484629810263795113e-77,
    9.6431887683992238428e-79,
    -4.2827980483017479213e-80,
    9.5087142369030441861e-82,
    2.7131392138694383464e-83,
    -4.0968779415069156659e-84,
    2.3742980019740160598e-85,
    -8.2770890210072789764e-87,
    9.072497609426645865e-89,
    1.0645558195026985633e-89,
    -9.285335619603754493e-91,
    4.3333135927203670323e-92,
    -1.1745606334673315984e-93,
    -2.6908010752365215433e-96,
    2.3898952892036810357e-96,
    -1.5569361182789167325e-97,
    6.0488748201074133757e-99,
    -1.2273370571029378615e-100,
    -2.540738850916238751e-102,
    3.7708800953170816508e-103,
    -2.0089261677502892352e-104,
    6.6158100911447349361e-106,
    -9.2404702022121568081e-108,
    -4.82072018655246532e-109,
    4.4938898756858357188e-110,
    -2.0497789059725778416e-111,
    5.7862770569866937508e-113,
    -4.5696744624334387424e-115,
    -5.8267365553303743945e-116,
    4.2025380699297338056e-117,
    -1.6889318527713702846e-118,
    4.1226213324018604871e-120,
    -8.2451196593745569675e-123,
    -5.2036993784470216679e-123,
    3.1616685922306712047e-124,
    -1.1432359131094236326e-125,
    2.4359648735131490197e-127,
    8.8701584767164321698e-130,
    -3.6328610892429035156e-130,
    1.9485148907440212068e-131,
    -6.450096583602651512e-133,
    1.215186561728963791e-134,
    1.0637863819629713691e-136,
    -2.0430980587447135517e-137,
    9.9760876002985183681e-139,
    -3.0707428945789381066e-140,
    5.2091832948433107534e-142,
    6.7131589510935005823e-144,
    -9.434301219575868381e-145,
    4.2908149482548296582e-146,
};

#define GAMMA_MIN_X 1.4616321449683623413
#define GAMMA_MIN_Y 0.88560319441088870028

/* Crude upper bound for psi(x) for x > 0, adequate for perturbation bounds
   for gamma. */
double
d_abs_digamma_ubound(double x)
{
    if (x <= 1.0)
    {
        return (1.0 + 1e-14) / x + 0.57721566490153286061 - x + 1e-14;
    }
    else if (x <= GAMMA_MIN_X)
    {
        return -1.250380137503405359*x + 1.8275958024049382196 + 1e-14;
    }
    else if (x <= 8.0)
    {
        return (x - GAMMA_MIN_X) * (1.7581621716802087234 +
            x * (-0.74622516195984912595 + x * (0.17009872711678924164 +
            x * (-0.018637559864260712285 + x * 0.00077747045691426195132)))) + 1e-12;
    }
    else if (x <= 128.0)
    {
        return 0.75334126757115431475 + x * (0.21045131598436795981 +
            x * (-0.0075387469533717503617 + x * (0.00017308475161765275722 +
                x * (-2.4025446500822043239e-6 + x * (1.9547402969088507111e-8 +
            x * (-8.5654894222045481692e-11 + x * 1.5584520745423393038e-13)))))) + 1e-12;
    }
    else
    {
        return (mag_d_log_upper_bound(x) + 1.0 / x) * (1.0 + 1e-14);
    }
}

/* Upper or lower bound (depending on direction) for gamma(x),
   assuming x > 0, no overflow. */
double
_arb_hypgeom_d_gamma(double x, int direction)
{
    double s, t, p;
    int i, r;

    if (direction == 1)
        p = 1 + 1e-14;
    else
        p = 1 - 1e-14;

    if (x < 0.5)
    {
        s = d_polyval(arb_hypgeom_rgamma_d_tab, 19, x);
        s = 1.0 / (s * x);
    }
    else if (x <= 1.5)
    {
        s = 1.0 / d_polyval(arb_hypgeom_rgamma_d_tab, 19, x - 1.0);
    }
    else
    {
        r = (int) (x + 0.5);

        s = d_polyval(arb_hypgeom_rgamma_d_tab, 19, x - r);

        t = 1.0;
        for (i = 0; i < r - 1; i++)
            t *= (x - i - 1) * p;

        s = t / s;
    }

    return s * p;
}

/* Set res = [a, b]; not checking overflow or underflow. */
void arb_set_interval_d_fast(arb_t res, double a, double b, slong prec)
{
    double mid, rad;

    if (a > b)
    {
        flint_throw(FLINT_ERROR, "arb_set_interval_d_fast: expected a < b\n");
    }

    mid = a + 0.5 * (b - a);
    rad = (0.5 * (b - a) + (mid * 1e-15)) * (1 + 1e-15);
    arf_set_d(arb_midref(res), mid);
    mag_set_d(arb_radref(res), rad);
    arb_set_round(res, res, prec);
}

int _arf_increment_fast(arf_t x, slong prec);



/* Try to compute gamma(x) using Taylor series. Returns 1 on success, 0 on
   failure (x too large or precision too large). */
int
arb_hypgeom_gamma_taylor(arb_t res, const arb_t x, int reciprocal, slong prec)
{
    double dx, dxerr, log2u, ds, du;
    slong i, n, wp, r, tail_bound, rad_exp, mid_exp;
    arf_t s, u, v;
    short term_prec[ARB_HYPGEOM_GAMMA_TAB_NUM];
    int success;

#if DEBUG
    printf("INPUT: "); arb_printd(x, 200); printf("\n");
    printf("INPUT prec: %ld\n", prec);
#endif

    /* We don't want to deal with infinities or huge/tiny exponents here. */
    if (!ARB_IS_LAGOM(x))
        return 0;

    /* 2^e bounds for the midpoint and radius. */
    mid_exp = arf_is_zero(arb_midref(x)) ? WORD_MIN : ARF_EXP(arb_midref(x));
    rad_exp = mag_is_zero(arb_radref(x)) ? WORD_MIN : MAG_EXP(arb_radref(x));

    /* Containing zero. */
    if (rad_exp >= mid_exp && arb_contains_zero(x))
    {
        if (reciprocal)
        {
            arb_t t;
            arb_init(t);
            arb_add_ui(t, x, 1, prec + 10);

            if (!arb_contains_zero(t))
            {
                success = arb_hypgeom_gamma_taylor(t, t, reciprocal, prec + 10);
                if (success)
                    arb_mul(res, x, t, prec);
            }
            else
            {
                /* todo: accurate wide interval */
                success = 0;
            }

            arb_clear(t);
            return success;
        }
        else
        {
            arb_indeterminate(res);
            return 1;
        }
    }

    /* Quick exclusion of too large numbers. */
    if (mid_exp > 8 || rad_exp > 8)
        return 0;

    /* Adjust precision if the input is not precise. */
    if (rad_exp != WORD_MIN)
        prec = FLINT_MIN(prec, -rad_exp + MAG_BITS);
    prec = FLINT_MAX(prec, 2);

    /* Midpoint and radius as doubles. */
    dx = arf_get_d(arb_midref(x), ARF_RND_NEAR);
    dxerr = mag_get_d(arb_radref(x));

    /* Too large to be efficient (high precision), or gamma(x) may overflow
       doubles (wide case). */
    if (dx + dxerr > 160.0 || dx - dxerr < -160.0)
        return 0;

    /* Very close to 0, reduce to gamma(x) = gamma(x + 1) / x. */
    if (mid_exp < -32 || (dx - dxerr >= -0.5 && dx - dxerr < ldexp(1.0, -6)))
    {
        arb_t t;
        arb_init(t);
        arb_add_ui(t, x, 1, prec + 10);

#if DEBUG
    printf("DIVIDING NEAR 0\n");
#endif

        success = arb_hypgeom_gamma_taylor(t, t, reciprocal, prec + 10);
        if (success)
        {
            if (reciprocal)
                arb_mul(res, x, t, prec);
            else
                arb_div(res, t, x, prec);
        }

        arb_clear(t);
        return success;
    }

    /* Nearest (roughly) integer to x, to use as shift for argument reduction
       to move to the interval [-0.5,0.5]. It's OK that dx is approximate so
       that the reduced argument will actually lie in [-0.5-eps,0.5+eps]. */
    if (dx >= 0.0)
        r = (slong) (dx + 0.5);
    else
        r = -(slong) (-dx + 0.5);

    /* Tuning cutoff. */
    if (prec >= 40)
    {
        if (r < -(40 + (prec - 40) / 4))
            return 0;

        if (r > 70 + (prec - 40) / 8)
            return 0;
    }

    /* For negative numbers, reduce to the positive case. */
    /* gamma(x) = (-1)^r * gamma(1+x-r) / (rf(1+r-x,-r)*(x-r)) */
    /* 1/gamma(x) = (-1)^r * rgamma(1+x-r) * rf(1+r-x,-r) * (x-r) */
    if (dx < 0.0)
    {
        arb_t t, u, v;

        arb_init(t);
        arb_init(u);
        arb_init(v);

        arb_sub_si(t, x, r, prec + 10);

        /* Pole. */
        if (!reciprocal && arb_contains_zero(t))
        {
            arb_indeterminate(res);
            success = 1;
        }
        else
        {
            arb_add_si(u, x, 1 - r, prec + 10);

            success = 1;
            if (reciprocal && !arb_is_positive(u))
            {
                /* todo: accurate wide interval */
                success = 0;
            }

            success = arb_hypgeom_gamma_taylor(u, u, reciprocal, prec + 10);

            if (success)
            {
                /* Wide bounds for rising factorial. */
                if (prec < 44)
                {
                    double a, b, c, d;

                    c = (-dx + r + 1 - dxerr) * (1 - 1e-14);
                    d = (-dx + r + 1 + dxerr) * (1 + 1e-14);
                    a = b = 1.0;

                    for (i = 0; i < -r; i++)
                    {
                        a = a * ((c + i) * (1 - 1e-15));
                        b = b * ((d + i) * (1 + 1e-15));
                    }

                    arb_set_interval_d_fast(v, a, b, 53);

                    if (reciprocal)
                    {
                        arb_mul(res, u, v, prec + 10);
                        arb_mul(res, res, t, prec);
                    }
                    else
                    {
                        arb_div(res, u, v, prec + 10);
                        arb_div(res, res, t, prec);
                    }
                }
                else
                {
                    arb_neg(v, x);
                    arb_add_si(v, v, 1 + r, prec + 10);
                    arb_hypgeom_rising_ui_rec(v, v, -r, prec + 10);
                    arb_mul(v, v, t, prec + 10);

                    if (reciprocal)
                        arb_mul(res, u, v, prec);
                    else
                        arb_div(res, u, v, prec);
                }

                if (r % 2)
                    arb_neg(res, res);
            }
        }

        arb_clear(t);
        arb_clear(u);
        arb_clear(v);
        return success;
    }

    /* Wide enclosure. */
    if (prec < 40 || rad_exp > -16)
    {
        double a, b, c;

#if DEBUG
    printf("WIDE CASE\n");
#endif

        dxerr += ldexp(1.0, mid_exp - 51);
        dxerr *= (1 + 1e-15);

        a = (dx - dxerr) * (1 - 1e-15);
        b = (dx + dxerr) * (1 + 1e-15);

        if (a >= GAMMA_MIN_X)
        {
            a = _arb_hypgeom_d_gamma(a, -1);
            b = _arb_hypgeom_d_gamma(b, 1);
        }
        else if (b <= GAMMA_MIN_X)
        {
            c = _arb_hypgeom_d_gamma(a, 1);
            a = _arb_hypgeom_d_gamma(b, -1);
            b = c;
        }
        else
        {
            a = _arb_hypgeom_d_gamma(a, 1);
            b = _arb_hypgeom_d_gamma(b, 1);
            b = FLINT_MAX(a, b);
            a = GAMMA_MIN_Y * (1 - 1e-15);
        }

        if (reciprocal)
        {
            c = (1.0 / b) * (1 - 1e-15);
            b = (1.0 / a) * (1 + 1e-15);
            a = c;
        }

        arb_set_interval_d_fast(res, a, b, prec);
        return 1;
    }

    /* Propagated error. */
    if (rad_exp == WORD_MIN)
    {
        dxerr = 0.0;
        rad_exp = WORD_MIN;
    }
    else
    {
        /* First-order relative error estimate plus safety factor to guarantee
           an upper bound. */
        dxerr = MAG_MAN(arb_radref(x)) * ldexp(1.0, -MAG_BITS);
        dxerr = dxerr * d_abs_digamma_ubound(dx) * 1.001;
    }

#if DEBUG
    flint_printf("propagated error = %g x 2^%wd\n", dxerr, rad_exp);
#endif

    wp = prec + 6 + FLINT_BIT_COUNT(FLINT_ABS(r));

    if (wp > ARB_HYPGEOM_GAMMA_TAB_PREC)
        return 0;

    success = 0;

    arf_init(s);
    arf_init(u);
    arf_init(v);

    /* u = x - r */
    arf_sub_si(u, arb_midref(x), r, wp, ARF_RND_DOWN);

    /* du = dx - r; */
    du = arf_get_d(u, ARF_RND_NEAR);

    /* bound log2(u) */
    if (-0.0001 < du && du < 0.0001)
        log2u = arf_is_zero(u) ? -wp : ARF_EXP(u);
    else
        log2u = mag_d_log_upper_bound(du < 0 ? -du : du) * 1.4426950408889634074 * (1 + 1e-14);

    term_prec[0] = wp;
    n = 0;

    for (i = 1; i < ARB_HYPGEOM_GAMMA_TAB_NUM; i++)
    {
        tail_bound = arb_hypgeom_gamma_coeffs[i].exp + i * log2u + 5;

        if (tail_bound <= -wp)
        {
            n = i;
            break;
        }

        term_prec[i] = FLINT_MIN(FLINT_MAX(wp + tail_bound, 2), wp);
    }

    if (n == 0)
    {
        flint_printf("warning: gamma_taylor: unexpected failure\n");
        success = 0;
        goto cleanup;
    }

#if DEBUG
    printf("COMPUTATION: wp = %ld, du = %g, log2u = %g, n = %ld\n", wp, du, log2u, n);
#endif

    if (wp <= 512 && n <= 128)
    {
        ds = 0.0;
        for (i = n - 1; i >= 1 && term_prec[i] <= 53; i--)
        {
#if DEBUG
            flint_printf("add term %wd with precision %wd (doubles)\n", i, term_prec[i]);
#endif

            ds = du * ds + arb_hypgeom_rgamma_d_tab[i];
        }

        arf_set_d(s, ds);
    }
    else
    {
        i = n - 1;
    }

    for ( ; i >= 1; i--)
    {
        arf_t c;

#if DEBUG
        flint_printf("add term %wd with precision %wd\n", i, term_prec[i]);
#endif

        if (!_arb_hypgeom_gamma_coeff_shallow(c, NULL, i, term_prec[i]))
        {
            flint_throw(FLINT_ERROR, "(arb_hypgeom_gamma_taylor):"
                    "prec = %wd, du = %g, log2u = %d, term_prec[%wd] = %wd",
                    prec, du, log2u, i, term_prec[i]);
        }

        if (term_prec[i] < wp - 128)
        {
            arf_set_round(v, u, term_prec[i], ARF_RND_DOWN);
            arf_mul(s, s, v, term_prec[i], ARF_RND_DOWN);
            arf_add(s, s, c, term_prec[i], ARF_RND_DOWN);
        }
        else
        {
            arf_mul(s, s, u, term_prec[i], ARF_RND_DOWN);
            arf_add(s, s, c, term_prec[i], ARF_RND_DOWN);
        }
    }

    if (i == 0)
    {
#if DEBUG
        flint_printf("add term %wd with precision %wd\n", i, term_prec[i]);
#endif

        arf_mul(s, s, u, wp, ARF_RND_DOWN);
        arf_add_ui(s, s, 1, wp, ARF_RND_DOWN);
    }

    if (r == 0 || r == 1)
    {
        if (r == 0)
            arf_mul(s, s, u, wp, ARF_RND_DOWN);

        if (reciprocal)
        {
            arf_set_round(arb_midref(res), s, prec, ARF_RND_DOWN);
        }
        else
        {
            arf_one(u);
            arf_div(arb_midref(res), u, s, prec, ARF_RND_DOWN);
        }
        arf_mag_set_ulp(arb_radref(res), arb_midref(res), prec - 1);
    }
    else if (wp <= 320 || r <= 3)
    {
        _arf_increment_fast(u, wp);
        arf_set(v, u);

        for (i = 2; i < r; i++)
        {
            _arf_increment_fast(u, wp);
            arf_mul(v, v, u, wp, ARF_RND_DOWN);
        }

        if (reciprocal)
            arf_div(arb_midref(res), s, v, prec, ARF_RND_DOWN);
        else
            arf_div(arb_midref(res), v, s, prec, ARF_RND_DOWN);
        arf_mag_set_ulp(arb_radref(res), arb_midref(res), prec - 1);
    }
    else
    {
        arb_t t;
        arb_init(t);
        _arf_increment_fast(u, wp);
        arb_set_arf(t, u);
        arb_hypgeom_rising_ui_rec(t, t, r - 1, wp);

        if (reciprocal)
        {
            arb_set_arf(res, s);
            arb_div(res, res, t, prec);
        }
        else
            arb_div_arf(res, t, s, prec);

        arf_mag_add_ulp(arb_radref(res), arb_radref(res), arb_midref(res), prec - 1);
        arb_clear(t);
    }

    /* Add propagated error. */
    if (dxerr != 0)
    {
        mag_t err;
        double dy;
        dy = arf_get_d(arb_midref(res), ARF_RND_UP);
        dxerr = dxerr * dy * (1 + 1e-15);
        MAG_SET_D_2EXP(MAG_MAN(err), MAG_EXP(err), dxerr, rad_exp);
        mag_add(arb_radref(res), arb_radref(res), err);
    }

    success = 1;

#if DEBUG
    printf("OUTPUT: "); arb_printd(res, 200); printf("\n");
#endif

cleanup:
    arf_clear(s);
    arf_clear(u);
    arf_clear(v);

    return success;
}

