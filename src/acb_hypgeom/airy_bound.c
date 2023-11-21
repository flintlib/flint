/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

void
arb_bound_exp_neg(mag_t b, const arb_t x)
{
    arb_t t;
    arb_init(t);
    arf_set_mag(arb_midref(t), arb_radref(x));
    arf_sub(arb_midref(t), arb_midref(x), arb_midref(t), MAG_BITS, ARF_RND_FLOOR);
    arf_neg(arb_midref(t), arb_midref(t));
    arb_exp(t, t, MAG_BITS);
    arb_get_mag(b, t);
    arb_clear(t);
}

/* Implements DLMF 9.7.17. We assume |zeta| >= 1/2 and |arg(z)| <= 2 pi/3 here,
   ignoring the smaller points which must be dealt with separately. */
void
acb_hypgeom_airy_bound_9_7_17(mag_t bound, const acb_t z, const acb_t zeta)
{
    mag_t D, t, u, v, zeta_lower, half;

    mag_init(D);
    mag_init(t);
    mag_init(u);
    mag_init(v);
    mag_init(zeta_lower);
    mag_init(half);

    mag_one(half);
    mag_mul_2exp_si(half, half, -1);

    acb_get_mag_lower(zeta_lower, zeta);
    mag_max(zeta_lower, zeta_lower, half);  /* by assumption */

    /* 2 chi(1) exp(7 pi / (72 |zeta|)) * c_1 */
    /* simplified bound assuming |zeta| >= 1/2 */
    mag_one(D);

    /* exp(-zeta) / (2 sqrt(pi)) * (1 + D / |zeta|) */

    /* t = exp(-zeta) / (2 sqrt(pi)) < exp(-zeta) * (73/256) */
    arb_bound_exp_neg(t, acb_realref(zeta));
    mag_mul_ui(t, t, 73);
    mag_mul_2exp_si(t, t, -8);

    /* u = 1 + D / |zeta| */
    mag_div(u, D, zeta_lower);
    mag_one(v);
    mag_add(u, u, v);

    mag_mul(bound, t, u);

    mag_clear(D);
    mag_clear(t);
    mag_clear(u);
    mag_clear(v);
    mag_clear(zeta_lower);
}

void
acb_hypgeom_airy_bound_arg_le_2pi3(mag_t A, mag_t B, const acb_t z, slong wp)
{
    acb_t zeta, z1;

    acb_init(zeta);
    acb_init(z1);

    acb_set_round(zeta, z, wp);
    acb_sqrt(zeta, zeta, wp);
    acb_cube(zeta, zeta, wp);
    acb_mul_2exp_si(zeta, zeta, 1);
    acb_div_ui(zeta, zeta, 3, wp);

    acb_hypgeom_airy_bound_9_7_17(A, z, zeta);

    /* Use Bi(z) = w_1 Ai(z) + 2 w_2 Ai(z exp(+/- 2pi i / 3)),
       where w_1, w_2 are roots of unity */
    if (B != NULL)
    {
        arb_sqrt_ui(acb_imagref(z1), 3, wp);
        arb_set_si(acb_realref(z1), -1);
        acb_mul_2exp_si(z1, z1, -1);

        /* multiply by exp(-2 pi i / 3) in upper half plane
           and by exp(2 pi i / 3) in lower half plane, to stay close
           to positive reals */

        /* in principle, we should be computing the union of both cases,
           but since Bi is conjugate symmetric with the maximum value
           attained on the positive real line, it's sufficient to consider
           the case centered on the midpoint */
        if (arf_sgn(arb_midref(acb_imagref(z))) >= 0)
            acb_conj(z1, z1);

        acb_mul(z1, z1, z, wp);
        acb_neg(zeta, zeta);   /* same effect regardless of exp(+/-2 pi i/3) */

        acb_hypgeom_airy_bound_9_7_17(B, z1, zeta);
        mag_mul_2exp_si(B, B, 1);
        mag_add(B, B, A);
    }

    acb_clear(zeta);
    acb_clear(z1);
}

/* near the negative half line, use
   |Ai(-z)|, |Bi(-z)| <= |Ai(z exp(pi i/3))| + |Ai(z exp(-pi i/3))| */
void
acb_hypgeom_airy_bound_arg_ge_2pi3(mag_t A, mag_t B, const acb_t z, slong wp)
{
    acb_t zeta, z1, z2;

    acb_init(zeta);
    acb_init(z1);
    acb_init(z2);

    arb_sqrt_ui(acb_imagref(z1), 3, wp);
    arb_one(acb_realref(z1));
    acb_mul_2exp_si(z1, z1, -1);
    acb_conj(z2, z1);

    acb_neg_round(zeta, z, wp);
    acb_mul(z1, z1, zeta, wp);

    acb_sqrt(zeta, zeta, wp);
    acb_cube(zeta, zeta, wp);
    acb_mul_2exp_si(zeta, zeta, 1);
    acb_div_ui(zeta, zeta, 3, wp);
    acb_mul_onei(zeta, zeta);

    acb_hypgeom_airy_bound_9_7_17(A, z1, zeta);

    /* conjugate symmetry */
    if (arb_is_zero(acb_imagref(z)))
    {
        mag_mul_2exp_si(A, A, 1);
    }
    else
    {
        mag_t D;
        mag_init(D);
        acb_mul(z2, z2, zeta, wp);
        acb_neg(zeta, zeta);
        acb_hypgeom_airy_bound_9_7_17(D, z2, zeta);
        mag_add(A, A, D);
        mag_clear(D);
    }

    if (B != NULL)
        mag_set(B, A);

    acb_clear(zeta);
    acb_clear(z1);
    acb_clear(z2);
}

static int
arg_lt_2pi3_fast(const acb_t z)
{
    arf_t t;
    mag_t x, y, s;
    int res;

    if (arb_is_zero(acb_imagref(z)) && arb_is_nonnegative(acb_realref(z)))
        return 1;

    arf_init(t);
    mag_init(x);
    mag_init(y);
    mag_init(s);

    arf_set_mag(t, arb_radref(acb_realref(z)));
    arf_sub(t, arb_midref(acb_realref(z)), t, MAG_BITS, ARF_RND_FLOOR);

    if (arf_sgn(t) >= 0)
    {
        res = 1;
    }
    else
    {
        arf_get_mag(x, t);
        arb_get_mag_lower(y, acb_imagref(z));
        mag_set_ui(s, 3);
        mag_sqrt(s, s);
        mag_mul(s, s, x);
        res = mag_cmp(s, y) <= 0;
    }

    arf_clear(t);
    mag_clear(x);
    mag_clear(y);
    mag_clear(s);

    return res;
}

static int
arg_gt_2pi3_fast(const acb_t z)
{
    arf_t t;
    mag_t x, y, s;
    int res;

    if (arb_is_zero(acb_imagref(z)) && arb_is_negative(acb_realref(z)))
        return 1;

    arf_init(t);
    mag_init(x);
    mag_init(y);
    mag_init(s);

    arf_set_mag(t, arb_radref(acb_realref(z)));
    arf_add(t, arb_midref(acb_realref(z)), t, MAG_BITS, ARF_RND_CEIL);

    if (arf_sgn(t) >= 0)
    {
        res = 0;
    }
    else
    {
        arf_get_mag_lower(x, t);
        arb_get_mag(y, acb_imagref(z));
        mag_set_ui_lower(s, 3);
        mag_sqrt_lower(s, s);
        mag_mul_lower(s, s, x);
        res = mag_cmp(s, y) >= 0;
    }

    arf_clear(t);
    mag_clear(x);
    mag_clear(y);
    mag_clear(s);

    return res;
}

void
acb_hypgeom_airy_bound(mag_t ai, mag_t aip, mag_t bi, mag_t bip, const acb_t z)
{
    acb_t zeta;
    mag_t A, B, D, zlo, zhi;
    slong wp;
    int near_zero;

    /* Fast and slightly tighter bounds for real z <= 0 */
    /* Todo: could implement bounds for z >= 0 too */
    if (acb_is_real(z) && arb_is_nonpositive(acb_realref(z)))
    {
        mag_init(zlo);
        mag_init(zhi);
        mag_init(A);
        mag_init(B);
        mag_init(D);

        if (ai != NULL || bi != NULL)
        {
            acb_get_mag_lower(zlo, z);
            mag_rsqrt(A, zlo);
            mag_sqrt(A, A);
            mag_mul_ui(A, A, 150);
            mag_set_ui(D, 160);
            mag_min(A, A, D);
            mag_mul_2exp_si(A, A, -8);

            if (ai != NULL) mag_set(ai, A);
            if (bi != NULL) mag_set(bi, A);
        }

        if (aip != NULL || bip != NULL)
        {
            acb_get_mag(zhi, z);
            mag_sqrt(A, zhi);
            mag_sqrt(A, A);
            mag_mul_ui(A, A, 150);
            mag_set_ui(D, 160);
            mag_max(B, A, D);
            mag_mul_2exp_si(B, B, -8);
            mag_set_ui(D, 67);
            mag_max(A, A, D);
            mag_mul_2exp_si(A, A, -8);

            if (aip != NULL) mag_set(aip, A);
            if (bip != NULL) mag_set(bip, B);
        }

        mag_clear(zlo);
        mag_clear(zhi);
        mag_clear(A);
        mag_clear(B);
        mag_clear(D);
        return;
    }

    acb_init(zeta);
    mag_init(A);
    mag_init(B);
    mag_init(D);
    mag_init(zlo);
    mag_init(zhi);

    wp = MAG_BITS * 2;

    acb_get_mag_lower(zlo, z);
    acb_get_mag(zhi, z);

    if (mag_cmp_2exp_si(zhi, 0) <= 0)
    {
        /* inside unit circle -- don't look at asymptotics */
        if (ai != NULL) mag_set_ui_2exp_si(ai, 159, -8);
        if (aip != NULL) mag_set_ui_2exp_si(aip, 125, -8);
        if (bi != NULL) mag_set_ui_2exp_si(bi, 310, -8);
        if (bip != NULL) mag_set_ui_2exp_si(bip, 239, -8);
    }
    else
    {
        /* look at asymptotics outside unit circle */
        near_zero = mag_cmp_2exp_si(zlo, 0) <= 0;
        if (near_zero)
            mag_one(zlo);

        if (arg_lt_2pi3_fast(z))
        {
            acb_hypgeom_airy_bound_arg_le_2pi3(A, (bi != NULL || bip != NULL) ? B : NULL, z, wp);
        }
        else if (arg_gt_2pi3_fast(z))
        {
            acb_hypgeom_airy_bound_arg_ge_2pi3(A, (bi != NULL || bip != NULL) ? B : NULL, z, wp);
        }
        else
        {
            mag_t A2, B2;

            mag_init(A2);
            mag_init(B2);

            acb_hypgeom_airy_bound_arg_le_2pi3(A, (bi != NULL || bip != NULL) ? B : NULL, z, wp);
            acb_hypgeom_airy_bound_arg_ge_2pi3(A2, (bi != NULL || bip != NULL) ? B2 : NULL, z, wp);

            mag_max(A, A, A2);
            mag_max(B, B, A2);

            mag_clear(A2);
            mag_clear(B2);
        }

        /* bound |z|^(1/4) */
        mag_sqrt(zhi, zhi);
        mag_sqrt(zhi, zhi);

        /* bound |z|^(-1/4) */
        mag_rsqrt(zlo, zlo);
        mag_sqrt(zlo, zlo);

        if (ai != NULL) mag_mul(ai, A, zlo);
        if (aip != NULL) mag_mul(aip, A, zhi);
        if (bi != NULL) mag_mul(bi, B, zlo);
        if (bip != NULL) mag_mul(bip, B, zhi);

        if (near_zero)
        {
            /* max by bounds on the unit circle */
            if (ai != NULL)
            {
                mag_set_ui_2exp_si(D, 159, -8);
                mag_max(ai, ai, D);
            }

            if (aip != NULL)
            {
                mag_set_ui_2exp_si(D, 125, -8);
                mag_max(aip, aip, D);
            }

            if (bi != NULL)
            {
                mag_set_ui_2exp_si(D, 310, -8);
                mag_max(bi, bi, D);
            }

            if (bip != NULL)
            {
                mag_set_ui_2exp_si(D, 239, -8);
                mag_max(bip, bip, D);
            }
        }
    }

    acb_clear(zeta);
    mag_clear(A);
    mag_clear(B);
    mag_clear(D);
    mag_clear(zlo);
    mag_clear(zhi);
}

