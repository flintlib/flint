/*
    Copyright (C) 2013, 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

/* r - |m| */
void
arb_get_mag_reverse(mag_t res, const arb_t x)
{
    mag_t t;
    mag_init(t);
    arf_get_mag_lower(t, arb_midref(x));
    mag_sub(res, arb_radref(x), t);
    mag_clear(t);
}

/* upper bound for re(rsqrt(x+yi)) / |rsqrt(x+yi)|,
   given upper bound for x, lower bound for y */
void
mag_rsqrt_re_quadrant1_upper(mag_t res, const mag_t x, const mag_t y)
{
    if (mag_is_zero(x))
    {
        mag_one(res);
        mag_mul_2exp_si(res, res, -1);
    }
    else
    {
        mag_t t, u;
        mag_init(t);
        mag_init(u);

        /* t = (y/x)^2 -- the result is a decreasing function of t */
        mag_div_lower(t, y, x);
        mag_mul_lower(t, t, t);

        /* (rsqrt(t^2+1)+1)/2 */
        mag_add_ui_lower(u, t, 1);
        mag_rsqrt(u, u);
        mag_add_ui(u, u, 1);
        mag_mul_2exp_si(res, u, -1);

        mag_clear(t);
        mag_clear(u);
    }

    mag_sqrt(res, res);
}

/* lower bound for re(rsqrt(x+yi)) / |rsqrt(x+yi)|,
   given lower bound for x, upper bound for y */
void
mag_rsqrt_re_quadrant1_lower(mag_t res, const mag_t x, const mag_t y)
{
    if (mag_is_zero(x))
    {
        mag_one(res);
        mag_mul_2exp_si(res, res, -1);
    }
    else
    {
        mag_t t, u;
        mag_init(t);
        mag_init(u);

        /* t = (y/x)^2 -- the result is a decreasing function of t */
        mag_div(t, y, x);
        mag_mul(t, t, t);

        /* (rsqrt(t^2+1)+1)/2 */
        mag_add_ui(u, t, 1);
        mag_rsqrt_lower(u, u);
        mag_add_ui_lower(u, u, 1);
        mag_mul_2exp_si(res, u, -1);

        mag_clear(t);
        mag_clear(u);
    }

    mag_sqrt_lower(res, res);
}

/* upper bound for re(rsqrt(-x+yi)) / |rsqrt(x+yi)|,
   given lower bound for -x, upper bound for y */
void
mag_rsqrt_re_quadrant2_upper(mag_t res, const mag_t x, const mag_t y)
{
    if (mag_is_zero(x))
    {
        mag_one(res);
        mag_mul_2exp_si(res, res, -1);
    }
    else
    {
        mag_t t, u, v;
        mag_init(t);
        mag_init(u);
        mag_init(v);

        /* t = (y/x)^2 -- the result is an increasing function of t */
        mag_div(t, y, x);
        mag_mul(t, t, t);

        /* t / (2*(t+1)*(rsqrt(t+1)+1)) */
        mag_add_ui(u, t, 1);
        mag_rsqrt_lower(v, u);
        mag_add_ui_lower(v, v, 1);
        mag_add_ui_lower(u, t, 1);
        mag_mul_lower(v, v, u);
        mag_mul_2exp_si(v, v, 1);
        mag_div(res, t, v);

        mag_clear(t);
        mag_clear(u);
        mag_clear(v);
    }

    mag_sqrt(res, res);
}

/* lower bound for re(rsqrt(-x+yi)) / |rsqrt(x+yi)|,
   given upper bound for -x, lower bound for y */
void
mag_rsqrt_re_quadrant2_lower(mag_t res, const mag_t x, const mag_t y)
{
    if (mag_is_zero(x))
    {
        mag_one(res);
        mag_mul_2exp_si(res, res, -1);
    }
    else
    {
        mag_t t, u, v;
        mag_init(t);
        mag_init(u);
        mag_init(v);

        /* t = (y/x)^2 -- the result is an increasing function of t */
        mag_div_lower(t, y, x);
        mag_mul_lower(t, t, t);

        /* t / (2*(t+1)*(rsqrt(t+1)+1)) */
        mag_add_ui_lower(u, t, 1);
        mag_rsqrt(v, u);
        mag_add_ui(v, v, 1);
        mag_add_ui(u, t, 1);
        mag_mul(v, v, u);
        mag_mul_2exp_si(v, v, 1);
        mag_div_lower(res, t, v);

        mag_clear(t);
        mag_clear(u);
        mag_clear(v);
    }

    mag_sqrt_lower(res, res);
}

void
acb_rsqrt_wide(acb_t res, const acb_t z, slong prec)
{
    mag_t ax, ay, bx, by, cx, cy, dx, dy, am, bm;
    mag_t one;

    mag_init(ax); mag_init(ay); mag_init(bx); mag_init(by);
    mag_init(cx); mag_init(cy); mag_init(dx); mag_init(dy);
    mag_init(am); mag_init(bm);
    mag_init(one);

    mag_one(one);

    /* magnitude */
    acb_get_mag(am, z);
    mag_rsqrt_lower(am, am);
    acb_get_mag_lower(bm, z);
    mag_rsqrt(bm, bm);

    /* upper or lower half plane */
    if (arb_is_nonnegative(acb_imagref(z)) || arb_is_negative(acb_imagref(z)))
    {
        if (arb_is_nonnegative(acb_realref(z)))
        {
            arb_get_mag_lower(ax, acb_realref(z));
            arb_get_mag(ay, acb_imagref(z));
            arb_get_mag(bx, acb_realref(z));
            arb_get_mag_lower(by, acb_imagref(z));

            mag_rsqrt_re_quadrant2_lower(cx, bx, by);
            mag_rsqrt_re_quadrant2_upper(dx, ax, ay);

            /* equivalent but more expensive than pythagoras
            mag_rsqrt_re_quadrant1_lower(ax, ax, ay);
            mag_rsqrt_re_quadrant1_upper(bx, bx, by);
            */

            mag_mul(ax, dx, dx);
            mag_sub_lower(ax, one, ax);
            mag_sqrt_lower(ax, ax);
            mag_mul_lower(bx, cx, cx);
            mag_sub(bx, one, bx);
            mag_sqrt(bx, bx);
        }
        else
        {
            if (arb_is_nonpositive(acb_realref(z)))
            {
                arb_get_mag(ax, acb_realref(z));
                arb_get_mag_lower(ay, acb_imagref(z));
                arb_get_mag_lower(bx, acb_realref(z));
                arb_get_mag(by, acb_imagref(z));

                /* equivalent but more expensive than pythagoras
                mag_rsqrt_re_quadrant1_lower(cx, bx, by);
                mag_rsqrt_re_quadrant1_upper(dx, ax, ay);
                */

                mag_rsqrt_re_quadrant2_lower(ax, ax, ay);
                mag_rsqrt_re_quadrant2_upper(bx, bx, by);
            }
            else if (arf_sgn(arb_midref(acb_realref(z))) >= 0)
            {
                arb_get_mag_reverse(ax, acb_realref(z));
                arb_get_mag_lower(ay, acb_imagref(z));
                arb_get_mag(bx, acb_realref(z));
                arb_get_mag_lower(by, acb_imagref(z));

                mag_rsqrt_re_quadrant2_lower(ax, ax, ay);
                mag_rsqrt_re_quadrant1_upper(bx, bx, by);
            }
            else
            {
                arb_get_mag(ax, acb_realref(z));
                arb_get_mag_lower(ay, acb_imagref(z));
                arb_get_mag_reverse(bx, acb_realref(z));
                arb_get_mag_lower(by, acb_imagref(z));

                mag_rsqrt_re_quadrant2_lower(ax, ax, ay);
                mag_rsqrt_re_quadrant1_upper(bx, bx, by);
            }

            /* pythagoras */
            mag_mul(cx, bx, bx);
            mag_sub_lower(cx, one, bx);
            mag_sqrt_lower(cx, cx);
            mag_mul_lower(dx, ax, ax);
            mag_sub(dx, one, dx);
            mag_sqrt(dx, dx);
        }

        mag_mul_lower(ax, ax, am);
        mag_mul_lower(cx, cx, am);
        mag_mul(bx, bx, bm);
        mag_mul(dx, dx, bm);

        if (arf_sgn(arb_midref(acb_imagref(z))) > 0)
        {
            arb_set_interval_mag(acb_realref(res), ax, bx, prec);
            arb_set_interval_mag(acb_imagref(res), cx, dx, prec);
            arf_neg(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(res)));
        }
        else
        {
            arb_set_interval_mag(acb_realref(res), ax, bx, prec);
            arb_set_interval_mag(acb_imagref(res), cx, dx, prec);
        }
    }
    else if (arb_is_positive(acb_realref(z)))
    {
        /* right half plane, straddling real line */
        int symmetric;

        symmetric = arf_is_zero(arb_midref(acb_imagref(z)));

        arb_get_mag_lower(ax, acb_realref(z));
        arb_get_mag(dy, acb_imagref(z));
        arb_get_mag_reverse(cy, acb_imagref(z));

        if (!symmetric)
            mag_rsqrt_re_quadrant2_lower(cx, ax, cy);
        mag_rsqrt_re_quadrant2_upper(dx, ax, dy);

        mag_one(bx);
        /* mag_rsqrt_re_quadrant1_lower(ax, ax, dy); */
        mag_mul(ax, dx, dx);
        mag_sub_lower(ax, one, ax);
        mag_sqrt_lower(ax, ax);

        mag_mul_lower(ax, ax, am);
        mag_mul(bx, bx, bm);
        mag_mul(cx, cx, bm);
        mag_mul(dx, dx, bm);

        if (symmetric)
            arb_set_interval_neg_pos_mag(acb_imagref(res), dx, dx, prec);
        else if (arf_sgn(arb_midref(acb_imagref(z))) > 0)
            arb_set_interval_neg_pos_mag(acb_imagref(res), dx, cx, prec);
        else
            arb_set_interval_neg_pos_mag(acb_imagref(res), cx, dx, prec);

        arb_set_interval_mag(acb_realref(res), ax, bx, prec);
    }
    else   /* left half plane, straddling branch cut */
    {
        mag_zero(ax);
        arb_get_mag_lower(bx, acb_realref(z));
        arb_get_mag(by, acb_imagref(z));
        mag_rsqrt_re_quadrant2_upper(bx, bx, by);

        mag_mul_lower(ax, ax, am);
        mag_mul(bx, bx, bm);
        arb_set_interval_mag(acb_realref(res), ax, bx, prec);

        /* cx, dx = 1,1 */
        arb_set_interval_neg_pos_mag(acb_imagref(res), bm, bm, prec);
    }

    mag_clear(ax); mag_clear(ay); mag_clear(bx); mag_clear(by);
    mag_clear(cx); mag_clear(cy); mag_clear(dx); mag_clear(dy);
    mag_clear(am); mag_clear(bm);
    mag_clear(one);
}

void
acb_rsqrt_precise(acb_t y, const acb_t x, slong prec)
{
#define a acb_realref(x)
#define b acb_imagref(x)
#define c acb_realref(y)
#define d acb_imagref(y)

    arb_t r, t, u, v;
    slong wp;

    /* based on the identity sqrt(z) = sqrt(r) (z+r) / |z+r|: */
    /* 1/sqrt(a+bi) = (1/v)((a+r) - b*i), r = |a+bi|, v = sqrt(r*(b^2+(a+r)^2)) */

    wp = prec + 6;

    arb_init(r);
    arb_init(t);
    arb_init(u);
    arb_init(v);

    /* u = b^2, r = |a+bi| */
    arb_mul(t, a, a, wp);
    arb_mul(u, b, b, wp);
    arb_add(r, t, u, wp);
    arb_sqrtpos(r, r, wp);

    /* t = a+r, v = r*(b^2+(a+r)^2) */
    arb_add(t, r, a, wp);
    arb_mul(v, t, t, wp);
    arb_add(v, v, u, wp);
    arb_mul(v, v, r, wp);

    /* v = 1/sqrt(v) */
    arb_rsqrt(v, v, wp);

    arb_mul(c, t, v, prec);
    arb_mul(d, b, v, prec);
    arb_neg(d, d);

    arb_clear(r);
    arb_clear(t);
    arb_clear(u);
    arb_clear(v);

#undef a
#undef b
#undef c
#undef d
}

void
acb_rsqrt(acb_t y, const acb_t x, slong prec)
{
    slong acc;

#define a acb_realref(x)
#define b acb_imagref(x)
#define c acb_realref(y)
#define d acb_imagref(y)

    if (acb_contains_zero(x))
    {
        acb_indeterminate(y);
        return;
    }

    if (arb_is_zero(b))
    {
        if (arb_is_nonnegative(a))
        {
            arb_rsqrt(c, a, prec);
            arb_zero(d);
            return;
        }
        else if (arb_is_nonpositive(a))
        {
            arb_neg(d, a);
            arb_rsqrt(d, d, prec);
            arb_neg(d, d);
            arb_zero(c);
            return;
        }
    }

    if (arb_is_zero(a))
    {
        if (arb_is_nonnegative(b))
        {
            arb_mul_2exp_si(c, b, 1);
            arb_rsqrt(c, c, prec);
            arb_neg(d, c);
            return;
        }
        else if (arb_is_nonpositive(b))
        {
            arb_mul_2exp_si(c, b, 1);
            arb_neg(c, c);
            arb_rsqrt(c, c, prec);
            arb_set(d, c);
            return;
        }
    }

    acc = acb_rel_accuracy_bits(x);

    if (acc < 25)
    {
        acb_rsqrt_wide(y, x, prec);
    }
    else
    {
        if (arb_is_positive(a))
        {
            acb_rsqrt_precise(y, x, prec);
        }
        else if (arb_is_nonnegative(b))
        {
            acb_neg(y, x);
            acb_rsqrt_precise(y, y, prec);
            acb_div_onei(y, y);
        }
        else if (arb_is_negative(b))
        {
            acb_neg(y, x);
            acb_rsqrt_precise(y, y, prec);
            acb_mul_onei(y, y);
        }
        else
        {
            acb_rsqrt_wide(y, x, prec);
        }
    }

#undef a
#undef b
#undef c
#undef d
}

void
acb_rsqrt_analytic(acb_t res, const acb_t z, int analytic, slong prec)
{
    if (analytic && arb_contains_zero(acb_imagref(z)) &&
        !arb_is_positive(acb_realref(z)))
    {
        acb_indeterminate(res);
    }
    else
    {
        acb_rsqrt(res, z, prec);
    }
}

