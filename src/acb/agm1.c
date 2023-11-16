/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_poly.h"

void mag_agm(mag_t res, const mag_t x, const mag_t y);

/* Checks that |arg(z)| <= 3 pi / 4 */
static int
acb_check_arg(const acb_t z)
{
    mag_t re, im;
    int res;

    if (!arb_contains_negative(acb_realref(z)))
        return 1;

    mag_init(re);
    mag_init(im);

    arb_get_mag(re, acb_realref(z));
    arb_get_mag_lower(im, acb_imagref(z));

    res = mag_cmp(re, im) < 0;

    mag_clear(re);
    mag_clear(im);

    return res;
}

static void
sqrtmul(acb_t c, const acb_t a, const acb_t b, slong prec)
{
    if (arb_is_positive(acb_realref(a)) &&
        arb_is_positive(acb_realref(b)))
    {
        acb_mul(c, a, b, prec);
        acb_sqrt(c, c, prec);
    }
    else if (arb_is_nonnegative(acb_imagref(a)) &&
             arb_is_nonnegative(acb_imagref(b)))
    {
        acb_mul(c, a, b, prec);
        acb_neg(c, c);
        acb_sqrt(c, c, prec);
        acb_mul_onei(c, c);
    }
    else if (arb_is_nonpositive(acb_imagref(a)) &&
             arb_is_nonpositive(acb_imagref(b)))
    {
        acb_mul(c, a, b, prec);
        acb_neg(c, c);
        acb_sqrt(c, c, prec);
        acb_mul_onei(c, c);
        acb_neg(c, c);
    }
    else
    {
        acb_t d;
        acb_init(d);
        acb_sqrt(c, a, prec);
        acb_sqrt(d, b, prec);
        acb_mul(c, c, d, prec);
        acb_clear(d);
    }
}

static void
acb_agm_close_taylor(acb_t res, acb_t z, acb_t z2,
    const acb_t aplusb, const acb_t aminusb,
    const mag_t err, slong prec)
{
    acb_div(z, aminusb, aplusb, prec);
    acb_sqr(z, z, prec);
    acb_sqr(z2, z, prec);

    acb_mul_si(res, z2, -469, prec);
    acb_addmul_si(res, z, -704, prec);
    acb_mul(res, res, z2, prec);
    acb_addmul_si(res, z2, -1280, prec);
    acb_mul_2exp_si(z, z, 12);
    acb_sub(res, res, z, prec);
    acb_add_ui(res, res, 16384, prec);
    acb_mul_2exp_si(res, res, -15);

    acb_add_error_mag(res, err);

    acb_mul(res, res, aplusb, prec);
}

static void
acb_agm1_around_zero(acb_t res, const acb_t z, slong prec)
{
    mag_t a, b;
    mag_init(a);
    mag_init(b);
    mag_one(a);
    acb_get_mag(b, z);
    mag_agm(a, a, b);
    acb_zero(res);
    acb_add_error_mag(res, a);
    mag_clear(a);
    mag_clear(b);
}

void
acb_agm1_basecase(acb_t res, const acb_t z, slong prec)
{
    acb_t a, b, t, u;
    mag_t err, err2;
    int isreal;

    isreal = acb_is_real(z) && arb_is_nonnegative(acb_realref(z));

    if (isreal)
    {
        acb_init(a);
        acb_one(a);
        arb_agm(acb_realref(res), acb_realref(a), acb_realref(z), prec);
        arb_zero(acb_imagref(res));
        acb_clear(a);
        return;
    }

    if (acb_is_zero(z))
    {
        acb_zero(res);
        return;
    }

    if (acb_is_one(z))
    {
        acb_one(res);
        return;
    }

    if (!acb_check_arg(z))
    {
        acb_agm1_around_zero(res, z, prec);
        return;
    }

    acb_init(a);
    acb_init(b);
    acb_init(t);
    acb_init(u);
    mag_init(err);
    mag_init(err2);

    acb_one(a);
    acb_set_round(b, z, prec);

    while (1)
    {
        acb_sub(u, a, b, prec);

        if (acb_contains_zero(u))
        {
            /* Dupont's thesis, p. 87: |M(z) - a_n| <= |a_n - b_n| */
            acb_set(res, a);
            acb_get_mag(err, u);
            acb_add_error_mag(res, err);
            break;
        }

        acb_add(t, a, b, prec);

        acb_get_mag(err, u);
        acb_get_mag_lower(err2, t);
        mag_div(err, err, err2);
        mag_geom_series(err, err, 10);
        mag_mul_2exp_si(err, err, -6);

        if (mag_cmp_2exp_si(err, -prec) < 0)
        {
            acb_agm_close_taylor(res, a, b, t, u, err, prec);
            break;
        }

        acb_mul_2exp_si(t, t, -1);

        sqrtmul(u, a, b, prec);

        acb_swap(t, a);
        acb_swap(u, b);
    }

    acb_clear(a);
    acb_clear(b);
    acb_clear(t);
    acb_clear(u);
    mag_clear(err);
    mag_clear(err2);
}

/*
    Computes (M(z), M'(z)) using a finite difference.
    Assumes z exact, |arg(z)| <= 3 pi / 4.
*/
void
acb_agm1_deriv_diff(acb_t Mz, acb_t Mzp, const acb_t z, slong prec)
{
    mag_t err, t, C;
    fmpz_t rexp, hexp;
    acb_t u, v;
    slong wp, qexp;
    int isreal;

    if (!acb_is_exact(z) || !acb_is_finite(z) ||
          acb_is_zero(z) || !acb_check_arg(z))
    {
        acb_indeterminate(Mz);
        acb_indeterminate(Mzp);
        return;
    }

    isreal = acb_is_real(z) && arb_is_nonnegative(acb_realref(z));

    /*
       |M^(k)(z) / k!| <= C * D^k  where
       C = max(1, |z| + r),
       D = 1/r,
       and 0 < r < |z|

        M(z+h) - M(z-h)
       |--------------- - M'(z)|  <=  D^3 h^2 / (1 - D h)
              2h

        M(z+h) + M(z-h)
       |--------------- - M(z)|   <=  D^2 h^2 / (1 - D h)
               2

       h D < 1.
    */

    fmpz_init(hexp);
    fmpz_init(rexp);
    mag_init(err);
    mag_init(t);
    mag_init(C);
    acb_init(u);
    acb_init(v);

    /* choose r = 2^rexp such that r < |z| */
    acb_get_mag_lower(t, z);
    fmpz_sub_ui(rexp, MAG_EXPREF(t), 2);

    /* Choose h = r/q = 2^hexp = 2^(rexp-qexp)
       with qexp = floor(prec/2) + 5

       D = 1/r = 2^-rexp

       f(z)  error <= C D^2 h^2 / (1-Dh)
       f'(z) error <= C D^3 h^2 / (1-Dh)

       1/(1-Dh) < 2, hence:

       f(z)  error <  2 C D^2 h^2 = C 2^(1-2*qexp)
       f'(z) error <  2 C D^3 h^2 = C 2^(1-rexp-2*qexp)
    */

    /* C = max(1, |z| + r) */
    acb_get_mag(C, z);
    mag_one(t);
    mag_mul_2exp_fmpz(t, t, rexp);
    mag_add(C, C, t);
    mag_one(t);
    mag_max(C, C, t);

    qexp = prec / 2 + 5;
    /*
    if (fmpz_sgn(rexp) < 0)
        qexp += fmpz_bits(rexp);
    */

    /* compute h = 2^hexp */
    fmpz_sub_ui(hexp, rexp, qexp);

    /* compute finite differences */
    wp = prec + qexp + 5;

    acb_one(u);
    acb_mul_2exp_fmpz(u, u, hexp);
    acb_add(u, z, u, wp);
    acb_agm1_basecase(u, u, wp);

    acb_one(v);
    acb_mul_2exp_fmpz(v, v, hexp);
    acb_sub(v, z, v, wp);
    acb_agm1_basecase(v, v, wp);

    acb_add(Mz, u, v, prec);
    acb_sub(Mzp, u, v, prec);
    acb_mul_2exp_si(Mz, Mz, -1);
    acb_mul_2exp_si(Mzp, Mzp, -1);
    fmpz_neg(hexp, hexp);
    acb_mul_2exp_fmpz(Mzp, Mzp, hexp);

    /* add error */
    mag_mul_2exp_si(err, C, 1 - 2 * qexp);

    if (isreal)
        arb_add_error_mag(acb_realref(Mz), err);
    else
        acb_add_error_mag(Mz, err);

    fmpz_neg(rexp, rexp);
    mag_mul_2exp_fmpz(err, err, rexp);

    if (isreal)
        arb_add_error_mag(acb_realref(Mzp), err);
    else
        acb_add_error_mag(Mzp, err);

    fmpz_clear(hexp);
    fmpz_clear(rexp);
    mag_clear(err);
    mag_clear(t);
    mag_clear(C);
    acb_clear(u);
    acb_clear(v);
}

/*
For input z + eps

First derivative bound: max(1, |z|+|eps|+r) / r
Second derivative bound: 2 max(1, |z|+|eps|+r) / r^2

This is assuming that the circle at z with radius |eps| + r
does not cross the negative half axis, which we check.
*/

void
acb_agm1_deriv_right(acb_t Mz, acb_t Mzp, const acb_t z, slong prec)
{
    if (acb_is_exact(z))
    {
        acb_agm1_deriv_diff(Mz, Mzp, z, prec);
    }
    else
    {
        if (!acb_is_finite(z) || !acb_check_arg(z))
        {
            acb_indeterminate(Mz);
            acb_indeterminate(Mzp);
        }
        else
        {
            acb_t t;
            mag_t r, eps, err, one;
            int isreal;

            acb_init(t);
            mag_init(r);
            mag_init(err);
            mag_init(one);
            mag_init(eps);

            isreal = acb_is_real(z) && arb_is_nonnegative(acb_realref(z));

            mag_hypot(eps, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));

            /* choose r avoiding overlap with negative half axis */
            if (arf_sgn(arb_midref(acb_realref(z))) < 0)
                arb_get_mag_lower(r, acb_imagref(z));
            else
                acb_get_mag_lower(r, z);

            mag_mul_2exp_si(r, r, -1);

            if (mag_is_zero(r))
            {
                acb_indeterminate(Mz);
                acb_indeterminate(Mzp);
            }
            else
            {
                acb_set(t, z);
                mag_zero(arb_radref(acb_realref(t)));
                mag_zero(arb_radref(acb_imagref(t)));

                acb_get_mag(err, z);
                mag_add(err, err, r);
                mag_add(err, err, eps);
                mag_one(one);
                mag_max(err, err, one);
                mag_mul(err, err, eps);

                acb_agm1_deriv_diff(Mz, Mzp, t, prec);

                mag_div(err, err, r);

                if (isreal)
                    arb_add_error_mag(acb_realref(Mz), err);
                else
                    acb_add_error_mag(Mz, err);

                mag_div(err, err, r);
                mag_mul_2exp_si(err, err, 1);

                if (isreal)
                    arb_add_error_mag(acb_realref(Mzp), err);
                else
                    acb_add_error_mag(Mzp, err);
            }

            acb_clear(t);
            mag_clear(r);
            mag_clear(err);
            mag_clear(one);
            mag_clear(eps);
        }
    }
}

void
acb_agm1(acb_t res, const acb_t z, slong prec)
{
    if (acb_is_zero(z))
    {
        acb_zero(res);
    }
    else if (!acb_is_finite(z))
    {
        acb_indeterminate(res);
    }
    else if (acb_contains_zero(z))
    {
        acb_agm1_around_zero(res, z, prec);
    }
    else if (arf_sgn(arb_midref(acb_realref(z))) >= 0)
    {
        acb_agm1_basecase(res, z, prec);
    }
    else if (acb_equal_si(z, -1))
    {
        acb_zero(res);
    }
    else
    {
        /* use M(1,z) = M((z+1)/2, sqrt(z))
                      = (z+1)/2 * M(1, 2 sqrt(z) / (z+1))
                      = sqrt(z) * M(1, (z+1) / (2 sqrt(z)) */
        acb_t t;

        acb_init(t);
        acb_add_ui(t, z, 1, prec);
        acb_mul_2exp_si(t, t, -1);

        if (acb_contains_zero(t))
        {
            mag_t ra, rb;

            mag_init(ra);
            mag_init(rb);

            acb_get_mag(ra, t);
            acb_get_mag(rb, z);
            mag_sqrt(rb, rb);

            mag_agm(ra, ra, rb);

            acb_zero(res);
            acb_add_error_mag(res, ra);

            mag_clear(ra);
            mag_clear(rb);
        }
        else if (acb_rel_accuracy_bits(t) > acb_rel_accuracy_bits(z))
        {
            acb_sqrt(res, z, prec);
            acb_div(res, res, t, prec);
            acb_agm1_basecase(res, res, prec);
            acb_mul(res, res, t, prec);
        }
        else
        {
            acb_sqrt(res, z, prec);
            acb_div(t, t, res, prec);
            acb_agm1_basecase(t, t, prec);
            acb_mul(res, res, t, prec);
        }

        acb_clear(t);
    }
}

void
acb_agm1_deriv(acb_t Mz, acb_t Mzp, const acb_t z, slong prec)
{
    /*
       u = 2 sqrt(z) / (1+z)

       Mz = (1+z) M(u) / 2
       Mzp = [M(u) - (z-1) M'(u) / ((1+z) sqrt(z))] / 2
    */

    if (arf_sgn(arb_midref(acb_realref(z))) >= 0)
    {
        if (acb_is_one(z))
        {
            acb_one(Mz);
            acb_mul_2exp_si(Mzp, Mz, -1);
        }
        else
            acb_agm1_deriv_right(Mz, Mzp, z, prec);
    }
    else
    {
        acb_t t, u, zp1, zm1;

        acb_init(t);
        acb_init(u);
        acb_init(zp1);
        acb_init(zm1);

        acb_sqrt(t, z, prec);
        acb_add_ui(zp1, z, 1, prec);
        acb_sub_ui(zm1, z, 1, prec);

        acb_div(u, t, zp1, prec);
        acb_mul_2exp_si(u, u, 1);

        acb_agm1_deriv_right(Mz, Mzp, u, prec);

        acb_mul(Mzp, Mzp, zm1, prec);
        acb_mul(t, t, zp1, prec);
        acb_div(Mzp, Mzp, t, prec);
        acb_sub(Mzp, Mz, Mzp, prec);
        acb_mul_2exp_si(Mzp, Mzp, -1);

        acb_mul(Mz, Mz, zp1, prec);
        acb_mul_2exp_si(Mz, Mz, -1);

        acb_clear(t);
        acb_clear(u);
        acb_clear(zp1);
        acb_clear(zm1);
    }
}

void
acb_agm1_cpx(acb_ptr m, const acb_t z, slong len, slong prec)
{
    if (len < 1)
        return;

    if (len == 1)
    {
        acb_agm1(m, z, prec);
        return;
    }

    if (len == 2)
    {
        acb_agm1_deriv(m, m + 1, z, prec);
        return;
    }

    if (len >= 3)
    {
        acb_t t, u, v;
        acb_ptr w;
        slong k, n;

        acb_init(t);
        acb_init(u);
        acb_init(v);
        w = _acb_vec_init(len);

        acb_agm1_deriv(w, w + 1, z, prec);

        /* invert series */
        acb_inv(w, w, prec);
        acb_mul(t, w, w, prec);
        acb_mul(w + 1, w + 1, t, prec);
        acb_neg(w + 1, w + 1);

        if (acb_is_one(z))
        {
            for (k = 2; k < len; k++)
            {
                n = k - 2;

                acb_mul_ui(w + k, w + n + 0, (n+1)*(n+1), prec);
                acb_addmul_ui(w + k, w + n + 1, 7+3*n*(3+n), prec);
                acb_div_ui(w + k, w + k, 2*(n+2)*(n+2), prec);
                acb_neg(w + k, w + k);
            }
        }
        else
        {
            /* t = 3z^2 - 1 */
            /* u = -1 / (z^3 - z) */
            acb_mul(t, z, z, prec);
            acb_mul(u, t, z, prec);
            acb_mul_ui(t, t, 3, prec);
            acb_sub_ui(t, t, 1, prec);
            acb_sub(u, u, z, prec);
            acb_inv(u, u, prec);
            acb_neg(u, u);

            /* use differential equation for second derivative */
            acb_mul(w + 2, z, w + 0, prec);
            acb_addmul(w + 2, t, w + 1, prec);
            acb_mul(w + 2, w + 2, u, prec);
            acb_mul_2exp_si(w + 2, w + 2, -1);

            /* recurrence */
            for (k = 3; k < len; k++)
            {
                n = k - 3;
                acb_mul_ui(w + k, w + n + 0, (n+1)*(n+1), prec);
                acb_mul(v, w + n + 1, z, prec);
                acb_addmul_ui(w + k, v, 7+3*n*(3+n), prec);
                acb_mul(v, w + n + 2, t, prec);
                acb_addmul_ui(w + k, v, (n+2)*(n+2), prec);
                acb_mul(w + k, w + k, u, prec);
                acb_div_ui(w + k, w + k, (n+2)*(n+3), prec);
            }
        }

        /* invert series */
        _acb_poly_inv_series(m, w, len, len, prec);

        acb_clear(t);
        acb_clear(u);
        acb_clear(v);
        _acb_vec_clear(w, len);
    }
}

