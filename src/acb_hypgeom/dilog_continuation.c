/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

static void
bsplit(acb_ptr VA, const acb_t z, const acb_t z2,
        const acb_t a, const acb_t a1a, slong k, slong h, slong prec)
{
    if (h - k == 1)
    {
        acb_zero(VA + 0);
        acb_mul_ui(VA + 1, a1a, (k+1)*(k+2), prec);
        acb_mul_si(VA + 2, z2, -k*k, prec);
        acb_mul_ui(VA + 3, a, (k+1)*(2*k+1), prec);
        acb_sub_ui(VA + 3, VA + 3, (k+1)*(k+1), prec);
        acb_mul(VA + 3, VA + 3, z, prec);
        acb_neg(VA + 3, VA + 3);
        acb_set(VA + 4, VA + 1);
        acb_zero(VA + 5);
        acb_set(VA + 6, VA + 1);
    }
    else
    {
        slong m;
        acb_ptr VB;

        if (h <= k)
            flint_throw(FLINT_ERROR, "(%s)\n", __func__);

        m = k + (h - k) / 2;
        VB = _acb_vec_init(7);

        bsplit(VA, z, z2, a, a1a, k, m, prec);
        bsplit(VB, z, z2, a, a1a, m, h, prec);

        acb_mul(VA + 6, VA + 6, VB + 6, prec);
        acb_mul(VA + 4, VA + 4, VB + 6, prec);
        acb_addmul(VA + 4, VA + 0, VB + 4, prec);
        acb_addmul(VA + 4, VA + 2, VB + 5, prec);
        acb_mul(VA + 5, VA + 5, VB + 6, prec);
        acb_addmul(VA + 5, VA + 1, VB + 4, prec);
        acb_addmul(VA + 5, VA + 3, VB + 5, prec);
        acb_set(VB + 6, VA + 3);
        acb_mul(VA + 3, VA + 3, VB + 3, prec);
        acb_addmul(VA + 3, VA + 1, VB + 2, prec);
        acb_set(VB + 5, VA + 2);
        acb_mul(VA + 2, VA + 2, VB + 3, prec);
        acb_addmul(VA + 2, VA + 0, VB + 2, prec);
        acb_mul(VA + 1, VA + 1, VB + 0, prec);
        acb_addmul(VA + 1, VB + 6, VB + 1, prec);
        acb_mul(VA + 0, VA + 0, VB + 0, prec);
        acb_addmul(VA + 0, VB + 5, VB + 1, prec);

        _acb_vec_clear(VB, 7);
    }
}

/*
Some possible approaches to bounding the Taylor coefficients c_k at
the expansion point a:

1. Inspection of the symbolic derivatives gives the trivial bound

    |c_k| <= (1+|log(1-a)|) / min(|a|,|a-1|)^k

   which is good enough when not close to 0.

2. Using Cauchy's integral formula, some explicit computation gives
   |c_k| <= 4/|1-a|^k when |a| <= 1/2 or a = +/- i. The constant
   could certainly be improved.

3. For k >= 1, c_k = 2F1(k,k,k+1,a) / k^2. Can we use monotonicity to get
   good estimates here when a is complex? Note that 2F1(k,k,k,a) = (1-a)^-k.

*/

void
acb_hypgeom_dilog_continuation(acb_t res, const acb_t a, const acb_t z, slong prec)
{
    acb_t za, a1, a1a, za2, log1a;
    acb_ptr V;
    slong n;
    double tr;
    mag_t tm, err, am;
    int real;

    if (acb_is_zero(a))
    {
        acb_hypgeom_dilog_zero_taylor(res, z, prec);
        return;
    }

    if (acb_eq(a, z))
    {
        acb_zero(res);
        return;
    }

    acb_init(za);
    acb_init(a1);
    acb_init(a1a);
    acb_init(za2);
    acb_init(log1a);
    mag_init(tm);
    mag_init(err);
    mag_init(am);

    acb_sub(za, z, a, prec);       /* z-a */
    acb_sub_ui(a1, a, 1, prec);    /* a-1  */
    acb_mul(a1a, a1, a, prec);     /* (a-1)a */
    acb_mul(za2, za, za, prec);    /* (z-a)^2 */
    acb_neg(log1a, a1);
    acb_log(log1a, log1a, prec);   /* log(1-a) */

    acb_get_mag(am, a);
    if (mag_cmp_2exp_si(am, -1) <= 0 ||
        (acb_is_exact(a) && arb_is_zero(acb_realref(a)) &&
            arf_cmpabs_ui(arb_midref(acb_imagref(a)), 1) == 0))
    {
        acb_get_mag_lower(am, a1);
        acb_get_mag(tm, za);
        mag_div(tm, tm, am);  /* tm = ratio */
        mag_set_ui(am, 4);    /* am = prefactor */
    }
    else
    {
        acb_get_mag_lower(am, a);
        acb_get_mag_lower(tm, a1);
        mag_min(am, am, tm);
        acb_get_mag(tm, za);
        mag_div(tm, tm, am);  /* tm = ratio */
        acb_get_mag(am, log1a);
        mag_add_ui(am, am, 1);  /* am = prefactor */
    }

    tr = mag_get_d_log2_approx(tm);
    if (tr < -0.1)
    {
        arf_srcptr rr, ii;

        rr = arb_midref(acb_realref(z));
        ii = arb_midref(acb_imagref(z));
        if (arf_cmpabs(ii, rr) > 0)
            rr = ii;

        /* target relative accuracy near 0 */
        if (arf_cmpabs_2exp_si(rr, -2) < 0 && arf_cmpabs_2exp_si(rr, -prec) > 0)
            n = (prec - arf_abs_bound_lt_2exp_si(rr)) / (-tr) + 1;
        else
            n = prec / (-tr) + 1;

        n = FLINT_MAX(n, 2);
    }
    else
    {
        n = 2;
    }

    mag_geom_series(err, tm, n);
    mag_mul(err, err, am);

    real = acb_is_real(a) && acb_is_real(z) &&
        arb_is_negative(acb_realref(a1)) &&
        mag_is_finite(err);

    if (n < 10)
    {
        /* forward recurrence - faster for small n and/or low precision, but
           must be avoided for large n since complex intervals blow up */
        acb_t s, t, u, v;
        slong k;

        acb_init(s);
        acb_init(t);
        acb_init(u);
        acb_init(v);

        acb_div(u, log1a, a, prec);
        acb_neg(u, u);

        if (n >= 3)
        {
            acb_inv(v, a1, prec);
            acb_add(v, v, u, prec);
            acb_div(v, v, a, prec);
            acb_mul_2exp_si(v, v, -1);
            acb_neg(v, v);
            acb_mul(v, v, za2, prec);
        }

        acb_mul(u, u, za, prec);
        acb_add(s, u, v, prec);

        for (k = 3; k < n; k++)
        {
            acb_mul_ui(u, u, (k - 2) * (k - 2), prec);
            acb_mul(u, u, za2, prec);
            acb_mul_ui(t, a, (k - 1) * (2 * k - 3), prec);
            acb_sub_ui(t, t, (k - 1) * (k - 1), prec);
            acb_mul(t, t, v, prec);
            acb_addmul(u, t, za, prec);
            acb_mul_ui(t, a1a, (k - 1) * k, prec);
            acb_neg(t, t);
            acb_div(u, u, t, prec);
            acb_add(s, s, u, prec);
            acb_swap(v, u);
        }

        acb_set(res, s);

        acb_clear(s);
        acb_clear(t);
        acb_clear(u);
        acb_clear(v);
    }
    else
    {
        /* binary splitting */
        V = _acb_vec_init(7);

        bsplit(V, za, za2, a, a1a, 1, 1 + n, prec);

        acb_mul(V + 1, V + 4, log1a, prec);
        acb_neg(V + 1, V + 1);
        acb_mul(V + 2, V + 5, za2, prec);
        acb_mul_2exp_si(V + 2, V + 2, -1);
        acb_mul(V + 1, V + 1, za, prec);
        acb_div(V + 3, V + 2, a1, prec);
        acb_sub(V + 1, V + 1, V + 3, prec);
        acb_div(V + 0, log1a, a, prec);
        acb_addmul(V + 1, V + 2, V + 0, prec);
        acb_mul(V + 6, V + 6, a, prec);

        acb_div(V + 0, V + 1, V + 6, prec);
        acb_set(res, V + 0);

        _acb_vec_clear(V, 7);
    }

    if (real)
        arb_add_error_mag(acb_realref(res), err);
    else
        acb_add_error_mag(res, err);

    acb_clear(za);
    acb_clear(a1);
    acb_clear(a1a);
    acb_clear(za2);
    acb_clear(log1a);
    mag_clear(tm);
    mag_clear(err);
    mag_clear(am);
}

