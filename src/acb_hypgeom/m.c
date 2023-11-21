/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"

#ifdef __GNUC__
# define fabs __builtin_fabs
# define sqrt __builtin_sqrt
#else
# include <math.h>
#endif

void
acb_hypgeom_m_asymp(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)
{
    acb_t t, u, v, c;

    acb_init(t);
    acb_init(u);
    acb_init(v);
    acb_init(c);

    acb_sub(c, b, a, prec);
    acb_neg(v, z);

    acb_hypgeom_u_asymp(t, a, b, z, -1, prec);
    acb_hypgeom_u_asymp(u, c, b, v, -1, prec);

    /* gamma(b-a) */
    acb_rgamma(v, c, prec);
    acb_mul(t, t, v, prec);

    /* z^(a-b) */
    acb_neg(c, c);
    acb_pow(v, z, c, prec);
    acb_mul(u, u, v, prec);

    /* gamma(a) */
    acb_rgamma(v, a, prec);
    acb_mul(u, u, v, prec);

    /* exp(z) */
    acb_exp(v, z, prec);
    acb_mul(u, u, v, prec);

    /* (-z)^(-a) */
    acb_neg(c, a);
    acb_neg(v, z);
    acb_pow(v, v, c, prec);
    acb_mul(t, t, v, prec);

    acb_add(t, t, u, prec);

    if (!regularized)
    {
        acb_gamma(v, b, prec);

        if (acb_is_finite(v))
            acb_mul(t, t, v, prec);
        else
            acb_indeterminate(t);
    }

    if (acb_is_real(a) && acb_is_real(b) && acb_is_real(z))
    {
        arb_zero(acb_imagref(t));
    }

    acb_swap(res, t);

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
    acb_clear(c);
}

void
_acb_hypgeom_m_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z,
    int regularized, slong prec, slong gamma_prec, int kummer)
{
    if (regularized)
    {
        /* Remove singularity */
        if (acb_is_int(b) && arb_is_nonpositive(acb_realref(b)) &&
            arf_cmpabs_2exp_si(arb_midref(acb_realref(b)), 30) < 0)
        {
            acb_t c, d, t, u;
            slong n;

            n = arf_get_si(arb_midref(acb_realref(b)), ARF_RND_DOWN);

            acb_init(c);
            acb_init(d);
            acb_init(t);
            acb_init(u);

            acb_sub(c, a, b, prec);
            acb_add_ui(c, c, 1, prec);

            acb_neg(d, b);
            acb_add_ui(d, d, 2, prec);

            _acb_hypgeom_m_1f1(t, c, d, z, 0, prec, gamma_prec, kummer);

            acb_pow_ui(u, z, 1 - n, prec);
            acb_mul(t, t, u, prec);

            acb_rising_ui(u, a, 1 - n, prec);
            acb_mul(t, t, u, prec);

            arb_fac_ui(acb_realref(u), 1 - n, prec);
            acb_div_arb(res, t, acb_realref(u), prec);

            acb_clear(c);
            acb_clear(d);
            acb_clear(t);
            acb_clear(u);
        }
        else
        {
            acb_t t;
            acb_init(t);
            acb_rgamma(t, b, gamma_prec);
            _acb_hypgeom_m_1f1(res, a, b, z, 0, prec, gamma_prec, kummer);
            acb_mul(res, res, t, prec);
            acb_clear(t);
        }
        return;
    }

    /* Kummer's transformation */
    if (kummer)
    {
        acb_t u, v;
        acb_init(u);
        acb_init(v);

        acb_sub(u, b, a, prec);
        acb_neg(v, z);

        _acb_hypgeom_m_1f1(u, u, b, v, regularized, prec, gamma_prec, 0);
        acb_exp(v, z, prec);
        acb_mul(res, u, v, prec);

        acb_clear(u);
        acb_clear(v);
        return;
    }

    if (acb_is_one(a))
    {
        acb_hypgeom_pfq_direct(res, NULL, 0, b, 1, z, -1, prec);
    }
    else
    {
        acb_struct c[3];
        c[0] = *a;
        c[1] = *b;

        acb_init(c + 2);
        acb_one(c + 2);

        acb_hypgeom_pfq_direct(res, c, 1, c + 1, 2, z, -1, prec);

        acb_clear(c + 2);
    }
}

void
acb_hypgeom_m_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)
{
    if (arf_sgn(arb_midref(acb_realref(z))) >= 0
        || (acb_is_int(a) && arb_is_nonpositive(acb_realref(a))))
    {
        _acb_hypgeom_m_1f1(res, a, b, z, regularized, prec, prec, 0);
    }
    else
    {
        _acb_hypgeom_m_1f1(res, a, b, z, regularized, prec, prec, 1);
    }
}

static double
hypotmx(double x, double y)
{
    if (x > 0.0 && x > 1e6 * fabs(y))
        return y * y / (2.0 * x);
    else
        return sqrt(x * x + y * y) - x;
}

void
acb_hypgeom_m_choose(int * asymp, int * kummer, slong * wp,
    const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)
{
    double x, y, t, cancellation;
    double input_accuracy, direct_accuracy, asymp_accuracy;
    slong m = WORD_MAX;
    slong n = WORD_MAX;

    if (acb_is_int(a) &&
            arf_cmpabs_2exp_si(arb_midref(acb_realref(a)), 30) < 0)
    {
        m = arf_get_si(arb_midref(acb_realref(a)), ARF_RND_DOWN);
    }

    if (acb_is_int(b) &&
            arf_cmpabs_2exp_si(arb_midref(acb_realref(b)), 30) < 0)
    {
        n = arf_get_si(arb_midref(acb_realref(b)), ARF_RND_DOWN);
    }

    *asymp = 0;
    *kummer = 0;
    *wp = prec;

    /* The 1F1 series terminates. */
    /* TODO: for large m, estimate extra precision here. */
    if (m <= 0 && m < n && m > -10 * prec && (n > 0 || !regularized))
    {
        *asymp = 0;
        return;
    }

    /* The 1F1 series terminates with the Kummer transform. */
    /* TODO: for large m, estimate extra precision here. */
    if (m >= 1 && n >= 1 && m < 0.1 * prec && n < 0.1 * prec && n <= m)
    {
        *asymp = 0;
        *kummer = 1;
        return;
    }

    input_accuracy = acb_rel_one_accuracy_bits(z);
    t = acb_rel_one_accuracy_bits(a);
    input_accuracy = FLINT_MIN(input_accuracy, t);
    t = acb_rel_one_accuracy_bits(b);
    input_accuracy = FLINT_MIN(input_accuracy, t);
    input_accuracy = FLINT_MAX(input_accuracy, 0.0);

    /* From here we ignore the values of a, b. Taking them into account is
       a possible future improvement... */

    /* Tiny |z|. */
    if ((arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 2) < 0 &&
         arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 2) < 0))
    {
        *asymp = 0;
        *wp = FLINT_MAX(2, FLINT_MIN(input_accuracy + 20, prec));
        return;
    }

    /* Huge |z|. */
    if ((arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 64) > 0 ||
         arf_cmpabs_2exp_si(arb_midref(acb_imagref(z)), 64) > 0))
    {
        *asymp = 1;
        *wp = FLINT_MAX(2, FLINT_MIN(input_accuracy + 20, prec));
        return;
    }

    x = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_DOWN);
    y = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_DOWN);

    asymp_accuracy = sqrt(x * x + y * y) * 1.44269504088896 - 5.0;

    /* The Kummer transformation gives less cancellation with the 1F1 series. */
    if (x < 0.0)
    {
        *kummer = 1;
        x = -x;
    }

    if (asymp_accuracy >= prec)
    {
        *asymp = 1;
        *wp = FLINT_MAX(2, FLINT_MIN(input_accuracy + 20, prec));
        return;
    }

    cancellation = hypotmx(x, y) * 1.44269504088896;

    direct_accuracy = input_accuracy - cancellation;

    if (direct_accuracy > asymp_accuracy)
    {
        *asymp = 0;
        *wp = FLINT_MAX(2, FLINT_MIN(input_accuracy + 20, prec + cancellation));
    }
    else
    {
        *asymp = 1;
        *wp = FLINT_MAX(2, FLINT_MIN(input_accuracy + 20, prec));
    }
}

void
acb_hypgeom_m_nointegration(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)
{
    int asymp, kummer;
    slong wp;

    acb_hypgeom_m_choose(&asymp, &kummer, &wp, a, b, z, regularized, prec);

    if (asymp)
    {
        acb_hypgeom_m_asymp(res, a, b, z, regularized, wp);
    }
    else
    {
        _acb_hypgeom_m_1f1(res, a, b, z, regularized, wp, FLINT_MIN(wp, prec), kummer);
    }

    acb_set_round(res, res, prec);
}

void
acb_hypgeom_m(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)
{
    acb_t res2;
    slong acc, max, t;

    acb_init(res2);

    acb_hypgeom_m_nointegration(res2, a, b, z, regularized, prec);

    acc = acb_rel_accuracy_bits(res2);

    if (acc < 0.5 * prec)
    {
        max = prec;
        t = acb_rel_accuracy_bits(z);
        max = FLINT_MIN(max, t);
        t = acb_rel_accuracy_bits(a);
        max = FLINT_MIN(max, t);
        t = acb_rel_accuracy_bits(b);
        max = FLINT_MIN(max, t);

        if (max > 2 && acc < 0.5 * max)
        {
            if (acb_is_real(a) && acb_is_real(b) && acb_is_real(z) &&
                arf_cmpabs_2exp_si(arb_midref(acb_realref(a)), 60) < 0 &&
                arf_cmpabs_2exp_si(arb_midref(acb_realref(b)), 60) < 0 &&
                arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 60) < 0)
            {
                arb_hypgeom_1f1_integration(acb_realref(res),
                        acb_realref(a), acb_realref(b), acb_realref(z), regularized, prec);
                arb_zero(acb_imagref(res));

                if (acb_rel_accuracy_bits(res) > acb_rel_accuracy_bits(res2) ||
                    (acb_is_finite(res) && !acb_is_finite(res2)))
                {
                    acb_swap(res, res2);
                }
            }
        }
    }

    acb_swap(res, res2);
    acb_clear(res2);
}

void
acb_hypgeom_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)
{
    acb_hypgeom_m(res, a, b, z, regularized, prec);
}
