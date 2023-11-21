/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "arb_hypgeom.h"
#include "arb_fmpz_poly.h"

#ifdef __GNUC__
# define ldexp __builtin_ldexp
# define log __builtin_log
# define sqrt __builtin_sqrt
#else
# include <math.h>
#endif

/* Actually an enclosure for |x| < 0.99. */
static void
arb_erfinv_approx_tiny(arb_t res, const arb_t x, slong prec)
{
    arb_t t;
    mag_t err;
    arb_init(t);
    mag_init(err);
    arb_get_mag(err, x);
    mag_pow_ui(err, err, 3);
    arb_const_sqrt_pi(t, prec);
    arb_mul_2exp_si(t, t, -1);
    arb_mul(res, x, t, prec);
    arb_add_error_mag(res, err);
    arb_clear(t);
    mag_clear(err);
}

/* First terms of asymptotic expansion. */
/* http://www.ams.org/journals/mcom/1976-30-136/S0025-5718-1976-0421040-7/S0025-5718-1976-0421040-7.pdf */
static double erfinv_approx_big(double one_sub_x)
{
    double eta, l, y;

    eta = -log(one_sub_x*sqrt(3.1415926535897932));
    l = log(eta);
    y = sqrt(eta - 0.5*l + (0.25*l - 0.5)/eta + (1/16.*l*l - 3/8.*l + 7./8)/(eta*eta)
                 + (l*l*l/48 - 7*l*l/32 + 17*l/16 - 107./48)/(eta*eta*eta)
                + (l*l*l*l/128 - 23*l*l*l/192 + 29*l*l/32 - 31*l/8 + 1489/192.)/(eta*eta*eta*eta));
    return y;
}

/* First terms of asymptotic expansion, when a double can overflow. */
static void
arb_erfinv_approx_huge(arb_t res, const arb_t one_minus_x, slong prec)
{
    arb_t eta, l, y;
    fmpz c[5];
    arb_ptr poly;
    mag_t err;

    arb_init(eta);
    arb_init(l);
    arb_init(y);
    mag_init(err);
    poly = _arb_vec_init(5);

    arb_const_sqrt_pi(eta, prec);
    arb_mul(eta, eta, one_minus_x, prec);
    arb_log(eta, eta, prec);
    arb_neg(eta, eta);
    arb_log(l, eta, prec);
    arb_mul_ui(y, eta, 12, prec);
    arb_inv(y, y, prec);
    arb_mul_2exp_si(poly + 0, l, -1);
    arb_neg(poly + 0, poly + 0);
    c[0] = -2 * 3; c[1] = 3;
    _arb_fmpz_poly_evaluate_arb(poly + 1, c, 2, l, prec);
    c[0] = 14 * 9; c[1] = -6 * 9; c[2] = 9;
    _arb_fmpz_poly_evaluate_arb(poly + 2, c, 3, l, prec);
    c[0] = -214 * 18; c[1] = 102 * 18; c[2] = -21 * 18; c[3] = 2 * 18;
    _arb_fmpz_poly_evaluate_arb(poly + 3, c, 4, l, prec);
    c[0] = 2978 * 54; c[1] = -1488 * 54; c[2] = 348 * 54; c[3] = -46 * 54; c[4] = 3 * 54;
    _arb_fmpz_poly_evaluate_arb(poly + 4, c, 5, l, prec);
    _arb_poly_evaluate(res, poly, 5, y, prec);
    arb_add(res, res, eta, prec);
    arb_sqrt(res, res, prec);

    /* Should be an enclosure for 1-x < 1e-300 (but not proved). */
    arb_get_mag(err, res);
    mag_mul_2exp_si(err, err, -50);

    arb_clear(eta);
    arb_clear(l);
    arb_clear(y);
    mag_clear(err);
    _arb_vec_clear(poly, 5);
}

/* Adapted from https://people.maths.ox.ac.uk/gilesm/codes/erfinv/
   with permission */
/* Only good up to about 1 - 1e-15 */
/* Todo: use a good approximation for erfcinv */
static double erfinv_approx(double x, double one_sub_x)
{
    double w, p;

    w = -log(one_sub_x * (1.0 + x));

    if (w < 6.250000)
    {
        w = w - 3.125000;
        p =  -3.6444120640178196996e-21;
        p =   -1.685059138182016589e-19 + p*w;
        p =   1.2858480715256400167e-18 + p*w;
        p =    1.115787767802518096e-17 + p*w;
        p =   -1.333171662854620906e-16 + p*w;
        p =   2.0972767875968561637e-17 + p*w;
        p =   6.6376381343583238325e-15 + p*w;
        p =  -4.0545662729752068639e-14 + p*w;
        p =  -8.1519341976054721522e-14 + p*w;
        p =   2.6335093153082322977e-12 + p*w;
        p =  -1.2975133253453532498e-11 + p*w;
        p =  -5.4154120542946279317e-11 + p*w;
        p =    1.051212273321532285e-09 + p*w;
        p =  -4.1126339803469836976e-09 + p*w;
        p =  -2.9070369957882005086e-08 + p*w;
        p =   4.2347877827932403518e-07 + p*w;
        p =  -1.3654692000834678645e-06 + p*w;
        p =  -1.3882523362786468719e-05 + p*w;
        p =    0.0001867342080340571352 + p*w;
        p =  -0.00074070253416626697512 + p*w;
        p =   -0.0060336708714301490533 + p*w;
        p =      0.24015818242558961693 + p*w;
        p =       1.6536545626831027356 + p*w;
    }
    else if (w < 16.000000)
    {
        w = sqrt(w) - 3.250000;
        p =   2.2137376921775787049e-09;
        p =   9.0756561938885390979e-08 + p*w;
        p =  -2.7517406297064545428e-07 + p*w;
        p =   1.8239629214389227755e-08 + p*w;
        p =   1.5027403968909827627e-06 + p*w;
        p =   -4.013867526981545969e-06 + p*w;
        p =   2.9234449089955446044e-06 + p*w;
        p =   1.2475304481671778723e-05 + p*w;
        p =  -4.7318229009055733981e-05 + p*w;
        p =   6.8284851459573175448e-05 + p*w;
        p =   2.4031110387097893999e-05 + p*w;
        p =   -0.0003550375203628474796 + p*w;
        p =   0.00095328937973738049703 + p*w;
        p =   -0.0016882755560235047313 + p*w;
        p =    0.0024914420961078508066 + p*w;
        p =   -0.0037512085075692412107 + p*w;
        p =     0.005370914553590063617 + p*w;
        p =       1.0052589676941592334 + p*w;
        p =       3.0838856104922207635 + p*w;
    }
    else
    {
        w = sqrt(w) - 5.000000;
        p =  -2.7109920616438573243e-11;
        p =  -2.5556418169965252055e-10 + p*w;
        p =   1.5076572693500548083e-09 + p*w;
        p =  -3.7894654401267369937e-09 + p*w;
        p =   7.6157012080783393804e-09 + p*w;
        p =  -1.4960026627149240478e-08 + p*w;
        p =   2.9147953450901080826e-08 + p*w;
        p =  -6.7711997758452339498e-08 + p*w;
        p =   2.2900482228026654717e-07 + p*w;
        p =  -9.9298272942317002539e-07 + p*w;
        p =   4.5260625972231537039e-06 + p*w;
        p =  -1.9681778105531670567e-05 + p*w;
        p =   7.5995277030017761139e-05 + p*w;
        p =  -0.00021503011930044477347 + p*w;
        p =  -0.00013871931833623122026 + p*w;
        p =       1.0103004648645343977 + p*w;
        p =       4.8499064014085844221 + p*w;
    }

    return p*x;
}

/* floating-point approximation of erfinv(x), to be validated */
static void
arb_hypgeom_erfinv_guess(arb_t res, const arb_t x, const arb_t one_sub_x, slong extraprec)
{
    if (arf_cmpabs_2exp_si(arb_midref(x), -30) < 0)
    {
        arb_erfinv_approx_tiny(res, x, 128);
    }
    else if (arf_cmpabs_2exp_si(arb_midref(one_sub_x), -52) >= 0)
    {
        double y;

        y = erfinv_approx(arf_get_d(arb_midref(x), ARF_RND_NEAR),
                          arf_get_d(arb_midref(one_sub_x), ARF_RND_NEAR));

        arf_set_d(arb_midref(res), y);
        mag_set_d(arb_radref(res), ldexp(y, -50));
    }
    else if (arf_cmpabs_2exp_si(arb_midref(one_sub_x), -1000) >= 0)
    {
        double t, y;

        t = arf_get_d(arb_midref(one_sub_x), ARF_RND_NEAR);
        y = erfinv_approx_big(t);

        arf_set_d(arb_midref(res), y);
        mag_set_d(arb_radref(res), ldexp(y, -26 + 0.1 * log(t)));
    }
    else
    {
        arb_erfinv_approx_huge(res, one_sub_x, 30 + extraprec);
    }
}

void
arb_hypgeom_erfinv_precise(arb_t res, const arb_t x, const arb_t one_sub_x, int near_one, slong prec)
{
    slong wp;
    arb_t f, fprime, root, mid, t;
    slong extraprec, goal;
    int validated;

    if (arb_is_zero(x))
    {
        arb_zero(res);
        return;
    }

    arb_init(f);
    arb_init(fprime);
    arb_init(root);
    arb_init(mid);
    arb_init(t);

    goal = prec * 1.001 + 5;

    extraprec = fmpz_bits(ARF_EXPREF(arb_midref(one_sub_x)));
    extraprec += 15;

    /* Start with a guess */
    arb_hypgeom_erfinv_guess(root, x, one_sub_x, extraprec);
    validated = 0;

    while (!validated || arb_rel_accuracy_bits(root) <= goal)
    {
        /* We should get double the accuracy. */
        wp = arb_rel_accuracy_bits(root) * 2 + extraprec;
        /* But don't set the precision unrealistically high. */
        wp = FLINT_MIN(wp, 4 * (goal + extraprec));

        /* In case of quadratic convergence, avoid doing the
           penultimate iteration at higher precision than needed. */
        if (validated && wp < goal && wp > 0.7 * goal + 2 * extraprec)
            wp = goal / 2 + 2 * extraprec;

        arb_set(mid, root);
        mag_zero(arb_radref(mid));

        /* f(y) = erf(y) - x  OR  (1 - x) - erfc(y) */
        /* 1/f'(y) = exp(y^2) * sqrt(pi)/2 */
        if (near_one)
        {
            arb_hypgeom_erfc(f, mid, wp);
            arb_sub(f, one_sub_x, f, wp);
        }
        else
        {
            arb_hypgeom_erf(f, mid, wp);
            arb_sub(f, f, x, wp);
        }
        arb_sqr(fprime, root, wp);
        arb_exp(fprime, fprime, wp);
        arb_const_sqrt_pi(t, wp);
        arb_mul(fprime, fprime, t, wp);
        arb_mul_2exp_si(fprime, fprime, -1);
        arb_mul(t, f, fprime, wp);
        arb_sub(t, mid, t, wp);

        if (arb_contains_interior(root, t))
        {
            /* Interval Newton proves inclusion. */
            /* printf("newton %d -> 1  %ld: ", validated, wp); arb_printd(t, 50); printf("\n"); */
            validated = 1;
            arb_swap(root, t);
        }
        else
        {
            /* Try to improve the guess with a floating-point Newton step. */
            /* printf("newton %d -> 0  %ld: ", validated, wp); arb_printd(t, 50); printf("\n"); */
            arb_sqr(fprime, mid, wp);
            arb_exp(fprime, fprime, wp);
            arb_const_sqrt_pi(t, wp);
            arb_mul(fprime, fprime, t, wp);
            arb_mul_2exp_si(fprime, fprime, -1);
            arb_mul(t, f, fprime, wp);
            /* The Newton correction is a good guess for the error. */
            arb_get_mag(arb_radref(root), t);
            mag_mul_2exp_si(arb_radref(root), arb_radref(root), 1);
            arb_sub(t, mid, t, wp);
            arf_swap(arb_midref(root), arb_midref(t));
            /* Use more working precision to get out of any numerical difficulties */
            extraprec = extraprec * 1.05 + 10;
            validated = 0;
        }

        /* Something went wrong. Could fall back to bisection if we
           want to guarantee convergence? */
        if (extraprec > 10 * prec + 10000)
        {
            arb_indeterminate(root);
            break;
        }
    }

    arb_set_round(res, root, prec);

    arb_clear(f);
    arb_clear(fprime);
    arb_clear(root);
    arb_clear(mid);
    arb_clear(t);
}









static slong
arb_adjust_precision(slong prec, slong acc)
{
    acc = FLINT_MIN(acc, prec);
    acc = FLINT_MAX(acc, 0);
    prec = FLINT_MIN(prec, acc + MAG_BITS);
    prec = FLINT_MAX(prec, 2);
    return prec;
}



void
arb_hypgeom_erfinv(arb_t res, const arb_t x, slong prec)
{
    arb_t x1;
    int near_one;

    if (arb_is_zero(x))
    {
        arb_zero(res);
        return;
    }

    if (arf_sgn(arb_midref(x)) < 0)
    {
        arb_neg(res, x);
        arb_hypgeom_erfinv(res, res, prec);
        arb_neg(res, res);
        return;
    }

    if (arb_is_one(x))
    {
        arb_pos_inf(res);
        return;
    }

    arb_init(x1);
    near_one = ARF_EXP(arb_midref(x)) == 0;

    if (near_one)
    {
        arb_sub_ui(x1, x, 1, ARF_PREC_EXACT);
        arb_neg(x1, x1);
    }
    else
    {
        arb_sub_ui(x1, x, 1, prec + 30);
        arb_neg(x1, x1);
    }

    if (arb_is_positive(x1))
    {
        mag_t err;
        slong acc;
        arb_t xm;

        mag_init(err);
        arb_init(xm);

        /* Propagated error bound based on derivative. */
        /* erfinv'(x) <= (1/2) sqrt(pi) / (1 - |x|) */
        arb_get_mag_lower(err, x1);
        mag_inv(err, err);
        mag_mul(err, err, arb_radref(x));
        mag_mul_ui(err, err, 227);
        mag_mul_2exp_si(err, err, -8);

        acc = arb_rel_accuracy_bits(x);
        prec = arb_adjust_precision(prec, acc);

        arb_get_mid_arb(xm, x);

        if (near_one)
        {
            arb_sub_ui(x1, xm, 1, ARF_PREC_EXACT);
            arb_neg(x1, x1);
        }
        else
        {
            arb_sub_ui(x1, xm, 1, prec + 30);
            arb_neg(x1, x1);
        }

        arb_hypgeom_erfinv_precise(res, xm, x1, near_one, prec);
        arb_add_error_mag(res, err);

        mag_clear(err);
        arb_clear(xm);
    }
    else
    {
        arb_indeterminate(res);
    }

    arb_clear(x1);
}

void
arb_hypgeom_erfcinv(arb_t res, const arb_t x1, slong prec)
{
    arb_t x;

    if (arb_is_one(x1))
    {
        arb_zero(res);
        return;
    }

    arb_init(x);

    if (arf_cmp_d(arb_midref(x1), 0.01) <= 0 && arb_is_positive(x1))
    {
        mag_t err;
        slong acc;
        arb_t x1m, xm;

        mag_init(err);
        arb_init(x1m);
        arb_init(xm);

        /* Propagated error bound based on derivative. */
        /* erfinv'(x) <= (1/2) sqrt(pi) / (1 - |x|) */
        arb_get_mag_lower(err, x1);
        mag_inv(err, err);
        mag_mul(err, err, arb_radref(x1));
        mag_mul_ui(err, err, 227);
        mag_mul_2exp_si(err, err, -8);

        acc = arb_rel_accuracy_bits(x1);
        prec = arb_adjust_precision(prec, acc);

        arb_get_mid_arb(x1m, x1);
        arb_sub_ui(xm, x1m, 1, 2 * prec + 100);
        arb_neg(xm, xm);

        arb_hypgeom_erfinv_precise(res, xm, x1m, 1, prec);
        arb_add_error_mag(res, err);

        mag_clear(err);
        arb_clear(xm);
        arb_clear(x1m);
    }
    else
    {
        arb_sub_ui(x, x1, 1, 2 * prec + 100);
        arb_neg(x, x);
        arb_hypgeom_erfinv(res, x, prec);
    }

    arb_clear(x);
}
