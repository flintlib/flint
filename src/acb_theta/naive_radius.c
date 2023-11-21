/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"
#include "acb_theta.h"

/* Assuming a >= 0, return R2 such that x - (a/2)*log(x)\geq b for all
   x\geq R2, and R2 is close to the smallest possible */

static void
invert_lin_plus_log(arf_t R2, slong a, const arb_t b, slong prec)
{
    arb_t x, y, t;
    arf_t z;
    slong k;

    arb_init(x);
    arb_init(y);
    arb_init(t);
    arf_init(z);

    if (a == 0)
    {
        arb_get_ubound_arf(R2, b, prec);
        goto exit;
    }

    /* minimum is at x=a/2 */
    arb_set_si(x, a);
    arb_div_si(x, x, 2, prec);
    arb_log(y, x, prec);
    arb_mul(y, y, x, prec);
    arb_sub(y, x, y, prec);

    /* x = max(a, 2*(b - min) + a log 2) is always large enough; then iterate
       function a few times */
    arb_sub(y, b, y, prec);
    arb_const_log2(t, prec);
    arb_mul_2exp_si(t, t, -1);
    arb_mul_si(t, t, a, prec);
    arb_add(y, y, t, prec);
    arb_max(y, y, x, prec);
    arb_mul_si(x, y, 2, prec);
    arb_get_ubound_arf(z, x, prec);
    arb_set_arf(x, z);

    for (k = 0; k < 4; k++)
    {
        arb_log(y, x, prec);
        arb_mul_si(y, y, a, prec);
        arb_div_si(y, y, 2, prec);
        arb_add(x, b, y, prec);
        arb_get_ubound_arf(z, x, prec);
        arb_set_arf(x, z);
    }

    arb_get_ubound_arf(R2, x, prec);
    goto exit;

  exit:
    {
        arb_clear(x);
        arb_clear(y);
        arb_clear(t);
        arf_clear(z);
    }
}

void
acb_theta_naive_radius(arf_t R2, arf_t eps, const arb_mat_t C, slong ord, slong prec)
{
    slong g = arb_mat_nrows(C);
    slong lp = ACB_THETA_LOW_PREC;
    arb_t b, temp;
    arf_t cmp;
    slong k;

    arb_init(b);
    arb_init(temp);
    arf_init(cmp);

    arf_one(eps);
    arf_mul_2exp_si(eps, eps, -prec);
    arb_set_arf(b, eps);
    arb_mul_2exp_si(b, b, -2 * g - 2);

    /* Solve R2^((g-1)/2+ord) exp(-R2) \leq b */
    arb_log(b, b, lp);
    arb_neg(b, b);
    invert_lin_plus_log(R2, g - 1 + 2 * ord, b, lp);

    /* Max with 4, 2*ord for formula to be valid */
    arf_set_si(cmp, FLINT_MAX(4, 2 * ord));
    arf_max(R2, R2, cmp);

    /* Set error */
    arb_one(b);
    for (k = 0; k < g; k++)
    {
        arb_inv(temp, arb_mat_entry(C, k, k), lp);
        arb_add_si(temp, temp, 1, lp);
        arb_mul(b, b, temp, lp);
    }
    arb_mul_arf(b, b, eps, lp);
    arb_get_ubound_arf(eps, b, lp);

    arb_clear(b);
    arb_clear(temp);
    arf_clear(cmp);
}
