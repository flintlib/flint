/*
    Copyright (C) 2019 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_hypgeom.h"

/* need special case for integer l and z = 0 since the recurrence relations break down */
static void
_acb_hypgeom_coulomb_f_int_jet(acb_ptr F, const acb_t l, const acb_t eta, const acb_t z, slong len, slong prec)
{
    acb_poly_struct a[1];
    acb_poly_struct b[2];
    acb_poly_t zx, M, zxi;
    acb_t t;
    slong k;
    int real;

    if (arf_cmp_si(arb_midref(acb_realref(l)), -1) < 0)
    {
        _acb_vec_indeterminate(F, len);
        return;
    }

    /* https://fungrim.org/entry/2a2f18/ */
    /* F = C * (z+x)^(l+1) e^(-+ i (z+x)) M(l + 1 -+  i eta, 2l+2, +- 2 i (z+x)) */

    acb_poly_init(a);
    acb_poly_init(b);
    acb_poly_init(b + 1);
    acb_poly_init(zx);
    acb_poly_init(M);
    acb_poly_init(zxi);
    acb_init(t);

    acb_poly_set_coeff_acb(zx, 0, z);
    acb_poly_set_coeff_si(zx, 1, 1);

    acb_div_onei(t, eta);
    acb_add(t, t, l, prec);
    acb_add_ui(t, t, 1, prec);
    acb_poly_set_acb(a, t);

    acb_add_ui(t, l, 1, prec);
    acb_mul_2exp_si(t, t, 1);
    acb_poly_set_acb(b, t);

    acb_poly_one(b + 1);

    acb_onei(t);
    acb_mul_2exp_si(t, t, 1);
    acb_poly_scalar_mul(zxi, zx, t, prec);

    acb_hypgeom_pfq_series_direct(M, a, 1, b, 2, zxi, 1, -1, len, prec);

    acb_poly_scalar_mul_2exp_si(zxi, zxi, -1);
    acb_poly_neg(zxi, zxi);
    acb_poly_exp_series(zxi, zxi, len, prec);
    acb_poly_mullow(M, M, zxi, len, prec);

    acb_add_ui(t, l, 1, prec);
    acb_poly_pow_acb_series(zxi, zx, t, len, prec);

    acb_poly_mullow(M, M, zxi, len, prec);

    {
        /* C = 2^l exp((-pi eta + lu + lv)/2) */
        acb_t C, lu, lv;

        acb_init(C);
        acb_init(lu);
        acb_init(lv);

        acb_add_ui(lu, l, 1, prec);
        acb_mul_onei(t, eta);
        acb_add(lu, lu, t, prec);

        acb_add_ui(lv, l, 1, prec);
        acb_div_onei(t, eta);
        acb_add(lv, lv, t, prec);

        acb_lgamma(lu, lu, prec);
        acb_lgamma(lv, lv, prec);

        acb_const_pi(t, prec);
        acb_add(C, lu, lv, prec);
        acb_submul(C, t, eta, prec);
        acb_mul_2exp_si(C, C, -1);
        acb_exp(C, C, prec);

        acb_set_ui(t, 2);
        acb_pow(t, t, l, prec);
        acb_mul(C, C, t, prec);

        acb_poly_scalar_mul(M, M, C, prec);

        acb_clear(C);
        acb_clear(lu);
        acb_clear(lv);
    }

    real = acb_is_real(z) && acb_is_real(eta);

    for (k = 0; k < len; k++)
    {
        acb_poly_get_coeff_acb(F + k, M, k);
        if (real)
            arb_zero(acb_imagref(F + k));
    }

    acb_poly_clear(a);
    acb_poly_clear(b);
    acb_poly_clear(b + 1);
    acb_poly_clear(zx);
    acb_poly_clear(M);
    acb_poly_clear(zxi);
    acb_clear(t);
}

static void
_acb_hypgeom_coulomb_jet(acb_ptr F, acb_ptr G, acb_ptr Hpos, acb_ptr Hneg, const acb_t l, const acb_t eta, const acb_t z, slong len, slong prec)
{
    acb_t l1, t, R, S;

    if (len <= 0)
        return;

    if (len == 1)
    {
        acb_hypgeom_coulomb(F, G, Hpos, Hneg, l, eta, z, prec);
        return;
    }

    if (acb_contains_zero(z))
    {
        if (F != NULL)
        {
            if (acb_is_int(l))
                _acb_hypgeom_coulomb_f_int_jet(F, l, eta, z, len, prec);
            else
                _acb_vec_indeterminate(F, len);
        }

        if (G != NULL) _acb_vec_indeterminate(G, len);
        if (Hpos != NULL) _acb_vec_indeterminate(Hpos, len);
        if (Hneg != NULL) _acb_vec_indeterminate(Hneg, len);
        return;
    }


    acb_init(l1);
    acb_init(t);
    acb_init(R);
    acb_init(S);

    acb_add_ui(l1, l, 1, prec);

    acb_hypgeom_coulomb(F, G, Hpos, Hneg, l, eta, z, prec);

    /* todo: somehow recycle the gamma function values for the two evaluations? */
    acb_hypgeom_coulomb((F == NULL) ? NULL : (F + 1),
                        (G == NULL) ? NULL : (G + 1),
                        (Hpos == NULL) ? NULL : (Hpos + 1),
                        (Hneg == NULL) ? NULL : (Hneg + 1),
                        l1, eta, z, prec);

    /* First derivatives:
       https://fungrim.org/entry/a51a4b/, https://fungrim.org/entry/2fec14/ */

    /* R_l = (2l+1) C_l / C_{l+1} = sqrt(l+i eta) sqrt(l-i eta) / l */
    if (acb_is_real(l) && acb_is_real(eta) && arb_is_nonzero(acb_realref(eta)))
    {
        acb_mul(R, l1, l1, prec);
        acb_addmul(R, eta, eta, prec);
        acb_sqrt(R, R, prec);
    }
    else
    {
        acb_mul_onei(t, eta);
        acb_add(t, t, l1, prec);
        acb_sqrt(t, t, prec);
        acb_div_onei(R, eta);
        acb_add(R, R, l1, prec);
        acb_sqrt(R, R, prec);
        acb_mul(R, R, t, prec);
    }
    acb_div(R, R, l1, prec);

    acb_div(S, l1, z, prec);
    acb_div(t, eta, l1, prec);
    acb_add(S, S, t, prec);

    /* todo: fix regular F at origin */
    if (F != NULL)
    {
        acb_mul(F + 1, F + 1, R, prec);
        acb_neg(F + 1, F + 1);
        acb_addmul(F + 1, F, S, prec);
    }

    if (G != NULL)
    {
        acb_mul(G + 1, G + 1, R, prec);
        acb_neg(G + 1, G + 1);
        acb_addmul(G + 1, G, S, prec);
    }

    if (Hpos != NULL)
    {
        acb_mul(Hpos + 1, Hpos + 1, R, prec);
        acb_neg(Hpos + 1, Hpos + 1);
        acb_addmul(Hpos + 1, Hpos, S, prec);
    }

    if (Hneg != NULL)
    {
        acb_mul(Hneg + 1, Hneg + 1, R, prec);
        acb_neg(Hneg + 1, Hneg + 1);
        acb_addmul(Hneg + 1, Hneg, S, prec);
    }

    if (len >= 3)
    {
        acb_t q, q2, w, w2;

        acb_init(q);
        acb_init(q2);
        acb_init(w);
        acb_init(w2);

        acb_inv(w, z, prec);
        acb_mul(w2, w, w, prec);

        /* https://fungrim.org/entry/07a654/ */
        /* F''/2 = q F,   q = (2eta/z + l(l+1)/z^2 - 1)/2 */

        acb_mul(q, l, l1, prec);
        acb_mul(q, q, w2, prec);
        acb_mul_2exp_si(q2, eta, 1);
        acb_addmul(q, q2, w, prec);
        acb_sub_ui(q, q, 1, prec);
        acb_mul_2exp_si(q, q, -1);

        if (F != NULL) acb_mul(F + 2, F, q, prec);
        if (G != NULL) acb_mul(G + 2, G, q, prec);
        if (Hpos != NULL) acb_mul(Hpos + 2, Hpos, q, prec);
        if (Hneg != NULL) acb_mul(Hneg + 2, Hneg, q, prec);

        /* https://fungrim.org/entry/faa118/ */
        /* F'''/6 = (2qF' - q2 F)/6,   q2 = 2(eta + l(l+1)/z)/z^2 */
        if (len >= 4)
        {
            acb_mul_2exp_si(q, q, 1);
            acb_mul(q2, l, l1, prec);
            acb_mul(q2, q2, w, prec);
            acb_add(q2, q2, eta, prec);
            acb_mul_2exp_si(q2, q2, 1);
            acb_mul(q2, q2, w2, prec);

            if (F != NULL)
            {
                acb_mul(F + 3, F + 1, q, prec);
                acb_submul(F + 3, F + 0, q2, prec);
                acb_div_ui(F + 3, F + 3, 6, prec);
            }

            if (G != NULL)
            {
                acb_mul(G + 3, G + 1, q, prec);
                acb_submul(G + 3, G, q2, prec);
                acb_div_ui(G + 3, G + 3, 6, prec);
            }

            if (Hpos != NULL)
            {
                acb_mul(Hpos + 3, Hpos + 1, q, prec);
                acb_submul(Hpos + 3, Hpos, q2, prec);
                acb_div_ui(Hpos + 3, Hpos + 3, 6, prec);
            }

            if (Hneg != NULL)
            {
                acb_mul(Hneg + 3, Hneg + 1, q, prec);
                acb_submul(Hneg + 3, Hneg, q2, prec);
                acb_div_ui(Hneg + 3, Hneg + 3, 6, prec);
            }
        }

        /* https://fungrim.org/entry/eca10b/ */
        if (len >= 5)
        {
            slong k;
            acb_ptr s;
            s = _acb_vec_init(4);

            /*
            s4 = -(k(k+7) + 12)
            s3 = 2 (k(k+5) + 6) z
            s2 = (k(k+3) + z^2 - 2z*eta - l(l+1) + 2)
            s1 = 2(z-eta)
            s0 = 1
            F[k+4] = (s0*F[k] + s1*F[k+1] + s2*F[k+2] + s3*F[k+3]) / s4 / z^2
            */
            acb_sub(s + 1, z, eta, prec);
            acb_mul_2exp_si(s + 1, s + 1, 1);

            acb_mul(q, z, z, prec);
            acb_mul(q2, eta, z, prec);
            acb_mul_2exp_si(q2, q2, 1);
            acb_sub(q, q, q2, prec);
            acb_submul(q, l, l1, prec);

            for (k = 0; k + 4 < len; k++)
            {
                acb_mul_si(s + 3, z, 2 * (k * (k + 5) + 6), prec);
                acb_add_si(s + 2, q, k * (k + 3) + 2, prec);

                if (F != NULL)
                {
                    acb_dot(F + k + 4, F + k, 0, F + k + 1, 1, s + 1, 1, 3, prec);
                    acb_div_si(F + k + 4, F + k + 4, -(k * (k + 7) + 12), prec);
                    acb_mul(F + k + 4, F + k + 4, w2, prec);
                }

                if (G != NULL)
                {
                    acb_dot(G + k + 4, G + k, 0, G + k + 1, 1, s + 1, 1, 3, prec);
                    acb_div_si(G + k + 4, G + k + 4, -(k * (k + 7) + 12), prec);
                    acb_mul(G + k + 4, G + k + 4, w2, prec);
                }

                if (Hpos != NULL)
                {
                    acb_dot(Hpos + k + 4, Hpos + k, 0, Hpos + k + 1, 1, s + 1, 1, 3, prec);
                    acb_div_si(Hpos + k + 4, Hpos + k + 4, -(k * (k + 7) + 12), prec);
                    acb_mul(Hpos + k + 4, Hpos + k + 4, w2, prec);
                }

                if (Hneg != NULL)
                {
                    acb_dot(Hneg + k + 4, Hneg + k, 0, Hneg + k + 1, 1, s + 1, 1, 3, prec);
                    acb_div_si(Hneg + k + 4, Hneg + k + 4, -(k * (k + 7) + 12), prec);
                    acb_mul(Hneg + k + 4, Hneg + k + 4, w2, prec);
                }
            }

            _acb_vec_clear(s, 4);
        }

        acb_clear(q);
        acb_clear(q2);
        acb_clear(w);
        acb_clear(w2);
    }

    acb_clear(l1);
    acb_clear(t);
    acb_clear(R);
    acb_clear(S);
}

void
acb_hypgeom_coulomb_jet(acb_ptr F, acb_ptr G, acb_ptr Hpos, acb_ptr Hneg, const acb_t l, const acb_t eta, const acb_t z, slong len, slong prec)
{
    acb_t t;  /* to allow aliasing with z */
    acb_init(t);
    acb_set(t, z);
    _acb_hypgeom_coulomb_jet(F, G, Hpos, Hneg, l, eta, t, len, prec);
    acb_clear(t);
}

