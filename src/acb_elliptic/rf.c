/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_elliptic.h"

#ifdef __GNUC__
# define pow __builtin_pow
#else
# include <math.h>
#endif

static const unsigned short den_ratio_tab[512] = {
    1,1,10,7,12,11,26,1,136,19,2,23,20,1,58,31,
    16,1,74,1,164,43,2,47,56,1,106,1,4,59,122,1,
    32,67,2,71,292,1,2,79,24,83,2,1,356,1,2,1,
    1552,1,202,103,4,107,218,1,904,1,2,1,44,1,10,127,
    64,131,2,1,548,139,2,1,8,1,298,151,4,1,314,1,
    16,163,2,167,52,1,346,1,8,179,362,1,4,1,2,191,
    6176,1,394,199,4,1,2,1,8,211,2,1,4,1,2,223,
    16,227,458,1,932,1,2,239,1928,1,2,1,4,251,2,1,
    32896,1,2,263,4,1,538,271,8,1,554,1,1124,283,2,1,
    272,1,586,1,4,1,2,1,8,307,2,311,1252,1,634,1,
    32,1,2,1,4,331,2,1,2696,1,2,7,4,347,698,1,
    5648,1,2,359,76,1,2,367,8,1,746,1,4,379,2,383,
    64,1,778,1,4,1,794,1,3208,1,2,1,1636,1,2,1,
    16,419,842,1,4,1,2,431,3464,1,2,439,4,443,2,1,
    14368,1,2,1,1828,1,922,463,8,467,2,1,4,1,2,479,
    16,1,2,487,4,491,2,1,8,499,2,503,4,1,1018,1,
    256,1,2,1,2084,523,2,1,184,1,2,1,4,1,1082,1,
    16,547,2,1,4,1,1114,1,8,563,2,1,2276,571,2,1,
    18464,1,2,1,4,587,2,1,4744,1,2,599,2404,1,2,607,
    16,1,1226,1,2468,619,2,1,40,1,2,631,4,1,2,1,
    41024,643,2,647,4,1,1306,1,8,659,1322,1,4,1,2,1,
    10768,1,1354,1,4,683,2,1,8,691,2,1,4,1,1402,1,
    32,1,1418,1,4,1,2,719,8,1,2,727,12,1,1466,1,
    16,739,2,743,4,1,2,751,8,1,1514,1,3044,1,2,2,
    49216,1,1546,1,4,1,2,1,8,787,2,1,4,1,1594,1,
    16,1,2,1,3236,811,2,1,8,1,1642,823,4,827,1658,1,
    32,1,2,839,116,1,2,1,8,1,1706,1,3428,859,2,863,
    16,1,2,1,4,1,1754,1,7048,883,2,887,4,1,2,1,
    64,1,2,1,4,907,2,911,8,1,2,919,4,1,2,1,
    14864,1,2,1,3748,1,1882,1,8,947,2,1,3812,1,2,1,
    992,1,2,967,4,971,2,1,7816,1,2,983,4,1,2,991,
    16,1,1994,1,4,1,2,1,8072,1,2026,1,4,1019,2042,1,
};

void
acb_elliptic_rf_taylor_sum(acb_t res, const acb_t E2, const acb_t E3, slong nterms, slong prec)
{
    fmpz_t den, c, d, e;
    acb_ptr E2pow;
    arb_ptr E2powr;
    acb_t s;
    slong x, y, XMAX, YMAX, NMAX, N;
    int real;

    NMAX = nterms - 1;
    YMAX = NMAX / 3;
    XMAX = NMAX / 2;
    real = acb_is_real(E2) && acb_is_real(E3);

    fmpz_init(den);
    fmpz_init(c);
    fmpz_init(d);
    fmpz_init(e);
    acb_init(s);
    if (real)
    {
        E2powr = _arb_vec_init(XMAX + 1);
        E2pow = NULL;
        _arb_vec_set_powers(E2powr, acb_realref(E2), XMAX + 1, prec);
    }
    else
    {
        E2pow = _acb_vec_init(XMAX + 1);
        E2powr = NULL;
        _acb_vec_set_powers(E2pow, E2, XMAX + 1, prec);
    }

    /* Compute universal denominator. */
    fmpz_one(den);
    for (N = 1; N <= NMAX; N++)
        fmpz_mul_ui(den, den, den_ratio_tab[N]);

    /* Compute initial coefficient rf(1/2,y) / y! */
    fmpz_set(c, den);
    for (y = 0; y < YMAX; y++)
    {
        fmpz_mul_ui(c, c, 2 * y + 1);
        fmpz_divexact_ui(c, c, 2 * y + 2);
    }

    acb_zero(res);

    for (y = YMAX; y >= 0; y--)
    {
        acb_zero(s);

        if (y != YMAX)
        {
            fmpz_mul_ui(c, c, 2 * y + 2);
            fmpz_divexact_ui(c, c, 2 * y + 1);
        }

        fmpz_set(d, c);

        /* Use powers with respect to E2 */
        for (x = 0; x <= XMAX; x++)
        {
            N = 2 * x + 3 * y;

            if (N <= NMAX)
            {
                fmpz_divexact_ui(e, d, 2 * N + 1);

                if (x % 2 == 1)
                    fmpz_neg(e, e);

                if (x != 0 || y != 0)
                {
                    if (real)
                        arb_addmul_fmpz(acb_realref(s), E2powr + x, e, prec);
                    else
                        acb_addmul_fmpz(s, E2pow + x, e, prec);
                }

                fmpz_mul_ui(d, d, 2 * x + 2 * y + 1);
                fmpz_divexact_ui(d, d, 2 * x + 2);
            }
        }

        /* Horner with respect to E3 */
        acb_mul(res, res, E3, prec);
        acb_add(res, res, s, prec);
    }

    acb_div_fmpz(res, res, den, prec);
    acb_add_ui(res, res, 1, prec);

    fmpz_clear(den);
    fmpz_clear(c);
    fmpz_clear(d);
    fmpz_clear(e);
    acb_clear(s);
    if (real)
        _arb_vec_clear(E2powr, XMAX + 1);
    else
        _acb_vec_clear(E2pow, XMAX + 1);
}

void
acb_elliptic_rf(acb_t res, const acb_t x, const acb_t y, const acb_t z,
                    int flags, slong prec)
{
    acb_t xx, yy, zz, sx, sy, sz, t;
    acb_t X, Y, Z, E2, E3;
    mag_t err, err2, prev_err;
    slong k, wp, accx, accy, accz, order;

    if (!acb_is_finite(x) || !acb_is_finite(y) || !acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (acb_contains_zero(x) + acb_contains_zero(y) + acb_contains_zero(z) > 1)
    {
        acb_indeterminate(res);
        return;
    }

    acb_init(xx); acb_init(yy); acb_init(zz);
    acb_init(sx); acb_init(sy); acb_init(sz);
    acb_init(X); acb_init(Y); acb_init(Z); acb_init(E2); acb_init(E3);
    acb_init(t);
    mag_init(err);
    mag_init(err2);
    mag_init(prev_err);

    order = 5; /* will be set later */

    acb_set(xx, x);
    acb_set(yy, y);
    acb_set(zz, z);

    /* First guess precision based on the inputs. */
    /* This does not account for mixing. */
    accx = acb_rel_accuracy_bits(xx);
    accy = acb_rel_accuracy_bits(yy);
    accz = acb_rel_accuracy_bits(zz);
    accx = FLINT_MAX(accx, accy);
    accx = FLINT_MAX(accx, accz);
    if (accx < prec - 20)
        prec = FLINT_MAX(2, accx + 20);
    wp = prec + 10 + FLINT_BIT_COUNT(prec);

    /* Must do at least one iteration. */
    for (k = 0; k < prec; k++)
    {
        acb_sqrt(sx, xx, wp);
        acb_sqrt(sy, yy, wp);
        acb_sqrt(sz, zz, wp);

        acb_add(t, sy, sz, wp);
        acb_mul(t, t, sx, wp);
        acb_addmul(t, sy, sz, wp);

        acb_add(xx, xx, t, wp);
        acb_add(yy, yy, t, wp);
        acb_add(zz, zz, t, wp);

        acb_mul_2exp_si(xx, xx, -2);
        acb_mul_2exp_si(yy, yy, -2);
        acb_mul_2exp_si(zz, zz, -2);

        /* Improve precision estimate and set expansion order. */
        /* Should this done for other k also? */
        if (k == 0)
        {
            accx = acb_rel_accuracy_bits(xx);
            accy = acb_rel_accuracy_bits(yy);
            accz = acb_rel_accuracy_bits(zz);
            accx = FLINT_MAX(accx, accy);
            accx = FLINT_MAX(accx, accz);
            if (accx < prec - 20)
                prec = FLINT_MAX(2, accx + 20);
            wp = prec + 10 + FLINT_BIT_COUNT(prec);

            if (acb_is_real(xx) && acb_is_real(yy) && acb_is_real(zz))
                order = 2.1 * pow(prec, 0.4);
            else
                order = 2.5 * pow(prec, 0.4);
            order = FLINT_MIN(order, 500);
            order = FLINT_MAX(order, 2);
        }

        /* Close enough? Quick estimate based on |x-y|/|x| and |x-z|/|x| */
        /* We also terminate if there is no improvement. */
        acb_sub(t, xx, yy, wp);
        acb_get_mag(err, t);
        acb_sub(t, xx, zz, wp);
        acb_get_mag(err2, t);
        mag_max(err, err, err2);
        acb_get_mag_lower(err2, xx);
        mag_div(err, err, err2);
        mag_pow_ui(err, err, order);
        if (mag_cmp_2exp_si(err, -prec) < 0 ||
                (k > 2 && mag_cmp(err, prev_err) > 0))
            break;

        mag_set(prev_err, err);
    }

    /* X = 1-x/t, Y = 1-y/t, Z = -X-Y, t = (x+y+z)/3 */
    acb_add(t, xx, yy, wp);
    acb_add(t, t, zz, wp);
    acb_div_ui(t, t, 3, wp);

    acb_div(X, xx, t, wp);
    acb_sub_ui(X, X, 1, wp);
    acb_neg(X, X);

    acb_div(Y, yy, t, wp);
    acb_sub_ui(Y, Y, 1, wp);
    acb_neg(Y, Y);

    acb_add(Z, X, Y, wp);
    acb_neg(Z, Z);

    /* E2 = XY-Z^2, E3 = XYZ */
    acb_mul(E2, X, Y, wp);
    acb_mul(E3, E2, Z, wp);
    acb_submul(E2, Z, Z, wp);

    /*
    Crude bound for the coefficient of
    X^n1 Y^n2 Z^n3 with n1+n2+n3=n: 2*(9/8)^n.
    */

    /* Error bound. */
    acb_get_mag(err, X);
    acb_get_mag(err2, Y);
    mag_max(err, err, err2);
    acb_get_mag(err2, Z);
    mag_max(err, err, err2);

    mag_mul_ui(err, err, 9);
    mag_mul_2exp_si(err, err, -3);

    mag_geom_series(err, err, order);
    mag_mul_2exp_si(err, err, 1);

    acb_elliptic_rf_taylor_sum(sx, E2, E3, order, wp);

    if (acb_is_real(X) && acb_is_real(Y))
        arb_add_error_mag(acb_realref(sx), err);
    else
        acb_add_error_mag(sx, err);

    acb_rsqrt(t, t, wp);
    acb_mul(res, sx, t, prec);

    acb_clear(xx); acb_clear(yy); acb_clear(zz);
    acb_clear(sx); acb_clear(sy); acb_clear(sz);
    acb_clear(X); acb_clear(Y); acb_clear(Z); acb_clear(E2); acb_clear(E3);
    acb_clear(t);
    mag_clear(err);
    mag_clear(err2);
    mag_clear(prev_err);
}

