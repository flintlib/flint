/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpq.h"
#include "acb.h"
#include "acb_hypgeom.h"

void
arb_const_airy_ai0_eval(arb_t y, slong prec)
{
    arb_t t; fmpq_t v; arb_init(t); fmpq_init(v);
    arb_set_ui(y, 3);
    arb_root_ui(y, y, 3, prec + 5); arb_mul(y, y, y, prec + 5);
    fmpq_set_si(v, 2, 3); arb_gamma_fmpq(t, v, prec + 5);
    arb_mul(y, y, t, prec + 5);
    arb_inv(y, y, prec);
    arb_clear(t); fmpq_clear(v);
}

void
arb_const_airy_ai1_eval(arb_t y, slong prec)
{
    arb_t t; fmpq_t v; arb_init(t); fmpq_init(v);
    arb_set_ui(y, 3);
    arb_root_ui(y, y, 3, prec + 5);
    fmpq_set_si(v, 1, 3); arb_gamma_fmpq(t, v, prec + 5);
    arb_mul(y, y, t, prec + 5);
    arb_inv(y, y, prec); arb_neg(y, y);
    arb_clear(t); fmpq_clear(v);
}

void
arb_const_airy_bi0_eval(arb_t y, slong prec)
{
    arb_t t; fmpq_t v; arb_init(t); fmpq_init(v);
    arb_set_ui(y, 3);
    arb_root_ui(y, y, 6, prec + 5);
    fmpq_set_si(v, 2, 3); arb_gamma_fmpq(t, v, prec + 5);
    arb_mul(y, y, t, prec + 5);
    arb_inv(y, y, prec);
    arb_clear(t); fmpq_clear(v);
}

void
arb_const_airy_bi1_eval(arb_t y, slong prec)
{
    arb_t t; fmpq_t v; arb_init(t); fmpq_init(v);
    arb_set_ui(y, 3);
    arb_root_ui(y, y, 6, prec + 5);
    fmpq_set_si(v, 1, 3); arb_gamma_fmpq(t, v, prec + 5);
    arb_div(y, y, t, prec);
    arb_clear(t); fmpq_clear(v);
}

ARB_DEF_CACHED_CONSTANT(arb_const_airy_ai0, arb_const_airy_ai0_eval)
ARB_DEF_CACHED_CONSTANT(arb_const_airy_ai1, arb_const_airy_ai1_eval)
ARB_DEF_CACHED_CONSTANT(arb_const_airy_bi0, arb_const_airy_bi0_eval)
ARB_DEF_CACHED_CONSTANT(arb_const_airy_bi1, arb_const_airy_bi1_eval)

static void
acb_hypgeom_airy_0f1_sum_inner(acb_t s, acb_srcptr t, slong m, slong n, slong alpha, int real, slong prec)
{
    slong j, k;
    mp_limb_t c, chi, clo;

    acb_zero(s);

    /* not implemented (coefficient overflow) */
    if (FLINT_BITS == 32 && n > 37000)
    {
        acb_indeterminate(s);
        return;
    }

    c = 1;
    j = (n - 1) % m;

    for (k = n - 1; k >= 0; k--)
    {
        if (k != 0)
        {
            umul_ppmm(chi, clo, c, 3 * k + alpha);

            if (chi == 0)
                umul_ppmm(chi, clo, clo, k);

            if (chi != 0)
            {
                acb_div_ui(s, s, c, prec);
                c = 1;
            }
        }

        if (real)
            arb_addmul_ui(acb_realref(s), acb_realref(t + j), c, prec);
        else
            acb_addmul_ui(s, t + j, c, prec);

        if (k != 0)
        {
            c = c * k * (3 * k + alpha);

            if (j == 0)
            {
                acb_mul(s, s, t + m, prec);
                j = m - 1;
            }
            else
            {
                j--;
            }
        }
    }

    acb_div_ui(s, s, c, prec);
}

/* s1 = 0F1(1/3, z/3)
   s2 = 0F1(2/3, z/3)
   s4 = 0F1(4/3, z/3)
   s5 = 0F1(5/3, z/3) */
static void
acb_hypgeom_airy_0f1_sum(acb_t s1, acb_t s2, acb_t s4, acb_t s5, const acb_t z, slong n, int real, slong prec)
{
    acb_ptr t;
    slong m;

    m = 2 * n_sqrt(n);
    m = FLINT_MAX(m, 1);

    t = _acb_vec_init(m + 1);
    _acb_vec_set_powers(t, z, m + 1, prec);

    if (s1 != NULL) acb_hypgeom_airy_0f1_sum_inner(s1, t, m, n, -2, real, prec);
    if (s2 != NULL) acb_hypgeom_airy_0f1_sum_inner(s2, t, m, n, -1, real, prec);
    if (s4 != NULL) acb_hypgeom_airy_0f1_sum_inner(s4, t, m, n, +1, real, prec);
    if (s5 != NULL) acb_hypgeom_airy_0f1_sum_inner(s5, t, m, n, +2, real, prec);

    _acb_vec_clear(t, m + 1);
}

void
acb_hypgeom_airy_direct(acb_t ai, acb_t aip, acb_t bi, acb_t bip, const acb_t z, slong n, slong prec)
{
    mag_t err, wmag, tm;
    int is_real;
    acb_t s1, s2, s4, s5, t, u;
    arb_t ai0, ai1, bi0, bi1;

    mag_init(err);
    mag_init(wmag);
    mag_init(tm);

    acb_init(s1);
    acb_init(s2);
    acb_init(s4);
    acb_init(s5);
    acb_init(t);
    acb_init(u);

    arb_init(ai0);
    arb_init(ai1);
    arb_init(bi0);
    arb_init(bi1);

    n = FLINT_MAX(n, 2);
    is_real = acb_is_real(z);
    acb_get_mag(wmag, z);

    /*
    With w = z^3/9, the terms are bounded by 3 w^n / [(n-1)!]^2.

      3 w^n            w          w^2
    ---------- [ 1 +  ---  +   -------  + ....]
    ((n-1)!)^2        n^2      (n+1)^2
    */
    mag_pow_ui(wmag, wmag, 3);
    mag_div_ui(wmag, wmag, 9);
    mag_pow_ui(err, wmag, n);

    mag_div_ui(tm, wmag, n);
    mag_div_ui(tm, tm, n);
    mag_geom_series(tm, tm, 0);
    mag_mul(err, err, tm);

    mag_rfac_ui(tm, n - 1);
    mag_mul(tm, tm, tm);
    mag_mul(err, err, tm);
    mag_mul_ui(err, err, 3);

    acb_cube(t, z, prec);
    acb_div_ui(t, t, 3, prec);
    acb_hypgeom_airy_0f1_sum(
        (aip != NULL || bip != NULL) ? s1 : NULL,
        (ai != NULL || bi != NULL) ? s2 : NULL,
        (ai != NULL || bi != NULL) ? s4 : NULL,
        (aip != NULL || bip != NULL) ? s5 : NULL, t, n, is_real, prec);

    if (is_real)
    {
        arb_add_error_mag(acb_realref(s1), err);
        arb_add_error_mag(acb_realref(s2), err);
        arb_add_error_mag(acb_realref(s4), err);
        arb_add_error_mag(acb_realref(s5), err);
    }
    else
    {
        acb_add_error_mag(s1, err);
        acb_add_error_mag(s2, err);
        acb_add_error_mag(s4, err);
        acb_add_error_mag(s5, err);
    }

    if (ai != NULL || aip != NULL)
    {
        arb_const_airy_ai0(ai0, prec);
        arb_const_airy_ai1(ai1, prec);
    }

    if (bi != NULL || bip != NULL)
    {
        arb_const_airy_bi0(bi0, prec);
        arb_const_airy_bi1(bi1, prec);
    }

    /* support aliasing with z */
    acb_set(t, z);

    if (ai != NULL || bi != NULL)
    {
        acb_mul(u, s4, t, prec);

        if (ai != NULL)
        {
            acb_mul_arb(ai, s2, ai0, prec);
            acb_addmul_arb(ai, u, ai1, prec);
        }

        if (bi != NULL)
        {
            acb_mul_arb(bi, s2, bi0, prec);
            acb_addmul_arb(bi, u, bi1, prec);
        }
    }

    if (aip != NULL || bip != NULL)
    {
        acb_mul(u, t, t, prec);
        acb_mul_2exp_si(u, u, -1);
        acb_mul(u, u, s5, prec);

        if (aip != NULL)
        {
            acb_mul_arb(aip, s1, ai1, prec);
            acb_addmul_arb(aip, u, ai0, prec);
        }

        if (bip != NULL)
        {
            acb_mul_arb(bip, s1, bi1, prec);
            acb_addmul_arb(bip, u, bi0, prec);
        }
    }

    mag_clear(err);
    mag_clear(wmag);
    mag_clear(tm);

    acb_clear(s1);
    acb_clear(s2);
    acb_clear(s4);
    acb_clear(s5);
    acb_clear(t);
    acb_clear(u);

    arb_clear(ai0);
    arb_clear(ai1);
    arb_clear(bi0);
    arb_clear(bi1);
}

