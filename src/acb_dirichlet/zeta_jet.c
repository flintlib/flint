/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpq.h"
#include "acb_poly.h"
#include "acb_dirichlet.h"

static void
_acb_dirichlet_zeta_jet(acb_t t, const acb_t h, int deflate, slong len, slong prec)
{
    acb_t a;
    acb_init(a);
    acb_one(a);

    /* Use the reflection formula when Re(s) <= -4 (-2.5 at low precision,
       where the guard bits below are relatively costly). For larger
       Re(s) < 0, Euler-Maclaurin summation directly at s is faster (by
       5-60%, most at s = -1, -2, ... for len = 2) and loses fewer bits
       than the product of four series; it only needs some guard bits
       growing with -Re(s). */
    if (arf_cmp_d(arb_midref(acb_realref(h)), (prec < 128) ? -2.5 : -4.0) <= 0)
    {
        /* zeta(s) = (2*pi)**s * sin(pi*s/2) / pi * gamma(1-s) * zeta(1-s) */
        acb_t pi, hcopy;
        acb_ptr f, s1, s2, s3, s4, u;
        slong i;

        acb_init(pi);
        acb_init(hcopy);
        f = _acb_vec_init(2);
        s1 = _acb_vec_init(len);
        s2 = _acb_vec_init(len);
        s3 = _acb_vec_init(len);
        s4 = _acb_vec_init(len);
        u = _acb_vec_init(len);
        acb_set(hcopy, h);

        acb_const_pi(pi, prec);

        /* s1 = (2*pi)**s */
        acb_mul_2exp_si(pi, pi, 1);
        _acb_poly_pow_cpx(s1, pi, h, len, prec);
        acb_mul_2exp_si(pi, pi, -1);

        /* s2 = sin(pi*s/2) / pi */
        acb_set(f, h);
        acb_one(f + 1);
        acb_mul_2exp_si(f, f, -1);
        acb_mul_2exp_si(f + 1, f + 1, -1);
        _acb_poly_sin_pi_series(s2, f, 2, len, prec);
        _acb_vec_scalar_div(s2, s2, len, pi, prec);

        /* s3 = gamma(1-s) */
        acb_sub_ui(f, hcopy, 1, prec);
        acb_neg(f, f);
        acb_set_si(f + 1, -1);
        _acb_poly_gamma_series(s3, f, 2, len, prec);

        /* s4 = zeta(1-s) */
        acb_sub_ui(f, hcopy, 1, prec);
        acb_neg(f, f);
        _acb_poly_zeta_cpx_series(s4, f, a, 0, len, prec);
        for (i = 1; i < len; i += 2)
            acb_neg(s4 + i, s4 + i);

        _acb_poly_mullow(u, s1, len, s2, len, len, prec);
        _acb_poly_mullow(s1, s3, len, s4, len, len, prec);
        _acb_poly_mullow(t, u, len, s1, len, len, prec);

        /* add 1/(1-(s+t)) = 1/(1-s) + t/(1-s)^2 + ... */
        if (deflate)
        {
            acb_sub_ui(u, hcopy, 1, prec);
            acb_neg(u, u);
            acb_inv(u, u, prec);
            for (i = 1; i < len; i++)
                acb_mul(u + i, u + i - 1, u, prec);
            _acb_vec_add(t, t, u, len, prec);
        }

        acb_clear(pi);
        acb_clear(hcopy);
        _acb_vec_clear(f, 2);
        _acb_vec_clear(s1, len);
        _acb_vec_clear(s2, len);
        _acb_vec_clear(s3, len);
        _acb_vec_clear(s4, len);
        _acb_vec_clear(u, len);
    }
    else
    {
        slong wp = prec;

        if (arf_sgn(arb_midref(acb_realref(h))) < 0)
            wp += (slong) (-arf_get_d(arb_midref(acb_realref(h)), ARF_RND_DOWN)
                    * FLINT_BIT_COUNT(prec) * 0.7) + 2;

        _acb_poly_zeta_cpx_series(t, h, a, deflate, len, wp);

        /* trivial zeros zeta(-2n) = 0 (the reflection formula gives
           them exactly) */
        if (!deflate && acb_is_exact(h) && acb_is_int(h) && arf_sgn(arb_midref(acb_realref(h))) < 0 &&
            arf_is_int_2exp_si(arb_midref(acb_realref(h)), 1))
            acb_zero(t);

        _acb_vec_set_round(t, t, len, prec);
    }

    acb_clear(a);
}

/* todo: should adjust precision to input accuracy */
void
acb_dirichlet_zeta_jet(acb_t res, const acb_t s, int deflate, slong len, slong prec)
{
    double cutoff;

    if (len == 1 && deflate == 0)
    {
        acb_zeta(res, s, prec);
        return;
    }

    if (deflate == 0 && (arb_contains_zero(acb_imagref(s))
                    && arb_contains_si(acb_realref(s), 1)))
    {
        _acb_vec_indeterminate(res, len);
        return;
    }

    /* the approximate functional equation for exact rational numbers
       of small height */
    if (acb_is_exact(s) && arb_is_zero(acb_imagref(s)))
    {
        fmpq_t x;
        int done = 0;

        fmpq_init(x);
        if (_acb_dirichlet_get_small_fmpq(x, s) &&
            _acb_dirichlet_l_jet_use_afe(1, x, deflate, len, prec))
        {
            acb_dirichlet_l_fmpq_afe_jet(res, x, NULL, NULL, deflate, len, prec);
            done = _acb_vec_is_finite(res, len);
        }
        fmpq_clear(x);

        if (done)
            return;
    }

    if (len > 2 || deflate != 0)
    {
        _acb_dirichlet_zeta_jet(res, s, deflate, len, prec);
    }
    else
    {
        cutoff = 24.0 * prec * sqrt(prec);

        if (arb_is_exact(acb_realref(s)) &&
            arf_cmp_2exp_si(arb_midref(acb_realref(s)), -1) == 0)
            cutoff *= 2.5;
        else
            cutoff *= 4.0;

        if (arf_cmpabs_d(arb_midref(acb_imagref(s)), cutoff) >= 0 &&
            arf_cmpabs_d(arb_midref(acb_realref(s)), 10 + prec * 0.1) <= 0)
        {
            acb_dirichlet_zeta_jet_rs(res, s, len, prec);
        }
        else
        {
            _acb_dirichlet_zeta_jet(res, s, deflate, len, prec);
        }
    }
}
