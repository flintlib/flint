/*
    Copyright (C) 2019 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

static void
acb_hypgeom_coulomb_is_real(int * C, int * F, int * G, const acb_t l1, const acb_t eta, const acb_t z)
{
    *C = *F = *G = 0;

    if (acb_is_real(l1) && acb_is_real(eta))
    {
        if (arb_is_positive(acb_realref(l1)) || arb_is_nonzero(acb_realref(eta)))
        {
            *C = 1;
        }

        if (acb_is_real(z))
        {
            if (arb_is_positive(acb_realref(z)))
            {
                *F = *G = 1;
            }

            if (acb_is_int(l1))
                *F = 1;
        }
    }
}

void
_acb_hypgeom_coulomb(acb_t F, acb_t G, acb_t Hpos, acb_t Hneg, const acb_t l, const acb_t eta, const acb_t z, int asymp, slong prec)
{
    acb_t u, v, lu, lv, z1, z2, m, h, T1, T2, U1, U2, H1, H2, C, theta;
    int C_real, F_real, G_real;
    int want_U1, want_U2, cut;

    acb_init(u); acb_init(v); acb_init(lu); acb_init(lv);
    acb_init(z1); acb_init(z2); acb_init(m); acb_init(h);
    acb_init(T1); acb_init(T2); acb_init(U1); acb_init(U2);
    acb_init(H1); acb_init(H2); acb_init(C); acb_init(theta);

    acb_indeterminate(U1);
    acb_indeterminate(U2);

    /* z1 = 2iz, z2 = -2iz,  */
    acb_mul_onei(z1, z);
    acb_mul_2exp_si(z1, z1, 1);
    acb_neg(z2, z1);

    if (asymp == -1)
        asymp = acb_hypgeom_u_use_asymp(z1, prec);

    /* Need the union of both sides of the branch cut for G, H+, H-. */
    if (arb_is_nonnegative(acb_imagref(z)) || arb_is_negative(acb_imagref(z)) || arb_is_positive(acb_realref(z)))
        cut = 0;
    else
        cut = 1;

    want_U1 = want_U2 = 0;

    if (asymp)
    {
        want_U1 = want_U2 = 1;
    }
    else
    {
        if (G != NULL || Hpos != NULL || Hneg != NULL)
        {
            if (arf_sgn(arb_midref(acb_imagref(z))) >= 0)
                want_U1 = 1;
            else
                want_U2 = 1;

            if (cut)
                want_U1 = want_U2 = 1;
        }
    }

    /* m = l+1 */
    acb_add_ui(m, l, 1, prec);

    acb_hypgeom_coulomb_is_real(&C_real, &F_real, &G_real, m, eta, z);

    /* u = 1+l+i eta, v = 1+l-i eta */
    acb_mul_onei(u, eta);
    acb_add(u, u, m, prec);
    acb_div_onei(v, eta);
    acb_add(v, v, m, prec);

    /* lu = lgamma(u), v = lgamma(v) */
    acb_lgamma(lu, u, prec);

    if (C_real)
        acb_conj(lv, lu);
    else
        acb_lgamma(lv, v, prec);

    /* m = 2l+2 */
    acb_mul_2exp_si(m, m, 1);

    if (asymp)
    {
        if (want_U1 && want_U2 && G_real)
        {
            acb_hypgeom_u_asymp(U1, u, m, z2, -1, prec);
            acb_conj(U2, U1);
        }
        else
        {
            if (want_U1) acb_hypgeom_u_asymp(U1, u, m, z2, -1, prec);
            if (want_U2) acb_hypgeom_u_asymp(U2, v, m, z1, -1, prec);
        }
    }
    else
    {
        if (want_U1 && want_U2 && G_real)
        {
            acb_hypgeom_u(U1, u, m, z2, prec);
            acb_pow(h, z2, u, prec);
            acb_mul(U1, U1, h, prec);
            acb_conj(U2, U1);
        }
        else
        {
            if (want_U1)
            {
                acb_hypgeom_u(U1, u, m, z2, prec);
                acb_pow(h, z2, u, prec);
                acb_mul(U1, U1, h, prec);
            }

            if (want_U2)
            {
                acb_hypgeom_u(U2, v, m, z1, prec);
                acb_pow(h, z1, v, prec);
                acb_mul(U2, U2, h, prec);
            }
        }
    }

    /* C = exp((-pi eta + lu + lv)/2) */
    acb_const_pi(C, prec);
    acb_mul(C, C, eta, prec);
    acb_neg(C, C);

    if (C_real)
    {
        acb_mul_2exp_si(T1, lu, 1);
        arb_zero(acb_imagref(T1));
        acb_add(C, C, T1, prec);
    }
    else
    {
        acb_add(C, C, lu, prec);
        acb_add(C, C, lv, prec);
    }

    acb_mul_2exp_si(C, C, -1);

    /* https://fungrim.org/entry/1976e1/ */
    if (asymp)
    {
        /* T1 = exp(-(-iz + lv + u log(z1)) U1 */
        /* T2 = exp(-(+iz + lu + v log(z2)) U2 */
        acb_log(T1, z1, prec);
        acb_mul(T1, T1, u, prec);
        acb_add(T1, T1, lv, prec);
        acb_mul_2exp_si(z1, z1, -1);
        acb_sub(T1, T1, z1, prec);
        acb_mul_2exp_si(z1, z1, 1);
        acb_neg(T1, T1);
        acb_exp(T1, T1, prec);
        acb_mul(T1, T1, U1, prec);

        if (F_real)
        {
            acb_mul_2exp_si(F, T1, 1);
            arb_zero(acb_imagref(F));
        }
        else
        {
            acb_log(T2, z2, prec);
            acb_mul(T2, T2, v, prec);
            acb_add(T2, T2, lu, prec);
            acb_mul_2exp_si(z2, z2, -1);
            acb_sub(T2, T2, z2, prec);
            acb_mul_2exp_si(z2, z2, 1);
            acb_neg(T2, T2);
            acb_exp(T2, T2, prec);
            acb_mul(T2, T2, U2, prec);

            /* F = (T1 + T2) z C */
            acb_add(F, T1, T2, prec);
        }
    }
    else
    {
        /* C *= exp(-iz) */
        acb_div_onei(F, z);
        acb_add(C, C, F, prec);
        /* https://fungrim.org/entry/2a2f18/ */
        acb_hypgeom_m(F, v, m, z1, 1, prec);
    }

    if (acb_contains_zero(z))
    {
        acb_exp(C, C, prec);
        /* (2z)^l without logarithm */
        acb_mul_2exp_si(h, z, 1);
        acb_pow(h, h, l, prec);
        acb_mul(C, C, h, prec);

        /* h = log(2z) */
        acb_indeterminate(h);
    }
    else
    {
        /* h = log(2z) */
        acb_mul_2exp_si(h, z, 1);
        acb_log(h, h, prec);

        acb_addmul(C, h, l, prec);
        acb_exp(C, C, prec);
    }

    acb_mul(F, F, C, prec);
    if (F_real)
        arb_zero(acb_imagref(F));
    acb_mul(F, F, z, prec);

    if (G != NULL || Hpos != NULL || Hneg != NULL)
    {
        /* theta = z - eta h - 0.5 l pi + (lu - lv) / (2i) */
        acb_sub(theta, lu, lv, prec);
        acb_div_onei(theta, theta);
        acb_mul_2exp_si(theta, theta, -1);
        acb_const_pi(H1, prec);
        acb_mul_2exp_si(H1, H1, -1);
        acb_submul(theta, H1, l, prec);
        acb_submul(theta, eta, h, prec);
        acb_add(theta, theta, z, prec);

        /* H1 = exp(+i theta) U1, H2 = exp(-i theta) U2 */
        acb_mul_onei(H1, theta);
        acb_exp_invexp(H1, H2, H1, prec);
        acb_mul(H1, H1, U1, prec);
        acb_mul(H2, H2, U2, prec);

        if (G != NULL)
        {
            /* https://fungrim.org/entry/e2efbf/ */
            if (asymp && arb_is_positive(acb_realref(z)))
            {
                if (G_real)
                {
                    if (arf_sgn(arb_midref(acb_imagref(z))) >= 0)
                        acb_set(G, H1);
                    else
                        acb_set(G, H2);
                    arb_zero(acb_imagref(G));
                }
                else
                {
                    acb_add(G, H1, H2, prec);
                    acb_mul_2exp_si(G, G, -1);
                }
            }
            else
            {
                /* https://fungrim.org/entry/8027e8/ */
                acb_div_onei(u, F);
                acb_add(u, H1, u, prec);
                /* https://fungrim.org/entry/69e5fb/ */
                acb_mul_onei(v, F);
                acb_add(v, H2, v, prec);

                if (cut)
                    acb_union(G, u, v, prec);
                else if (arf_sgn(arb_midref(acb_imagref(z))) >= 0)
                    acb_set(G, u);
                else
                    acb_set(G, v);

                if (G_real)
                    arb_zero(acb_imagref(G));
            }
        }

        if (Hpos != NULL)
        {
            /* https://fungrim.org/entry/bcdfc6/ */
            acb_set(u, H1);
            /* https://fungrim.org/entry/f0414a/ */
            acb_mul_onei(v, F);
            acb_mul_2exp_si(v, v, 1);
            acb_add(v, H2, v, prec);

            if (cut)
                acb_union(Hpos, u, v, prec);
            else if (arf_sgn(arb_midref(acb_imagref(z))) >= 0)
                acb_set(Hpos, u);
            else
                acb_set(Hpos, v);

            if (G_real)
                arb_set(acb_imagref(Hpos), acb_realref(F));
        }

        if (Hneg != NULL)
        {
            /* https://fungrim.org/entry/0cc301/ */
            acb_div_onei(u, F);
            acb_mul_2exp_si(u, u, 1);
            acb_add(u, H1, u, prec);
            /* https://fungrim.org/entry/781eae/ */
            acb_set(v, H2);

            if (cut)
                acb_union(Hneg, u, v, prec);
            else if (arf_sgn(arb_midref(acb_imagref(z))) >= 0)
                acb_set(Hneg, u);
            else
                acb_set(Hneg, v);

            if (G_real)
                arb_neg(acb_imagref(Hneg), acb_realref(F));
        }
    }

    acb_clear(u); acb_clear(v); acb_clear(lu); acb_clear(lv);
    acb_clear(z1); acb_clear(z2); acb_clear(m); acb_clear(h);
    acb_clear(T1); acb_clear(T2); acb_clear(U1); acb_clear(U2);
    acb_clear(H1); acb_clear(H2); acb_clear(C); acb_clear(theta);
}

void
acb_hypgeom_coulomb(acb_t F, acb_t G, acb_t Hpos, acb_t Hneg, const acb_t l, const acb_t eta, const acb_t z, slong prec)
{
    /* We always compute F. Also handle aliasing. */
    acb_t F2, l2, eta2, z2;

    acb_init(F2);
    acb_init(l2);
    acb_init(eta2);
    acb_init(z2);

    acb_set(l2, l);
    acb_set(eta2, eta);
    acb_set(z2, z);

    _acb_hypgeom_coulomb(F2, G, Hpos, Hneg, l2, eta2, z2, -1, prec);

    if (F != NULL)
        acb_set(F, F2);

    acb_clear(F2);
    acb_clear(l2);
    acb_clear(eta2);
    acb_clear(z2);
}

