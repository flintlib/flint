/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"

#ifdef __GNUC__
# define sqrt __builtin_sqrt
#else
# include <math.h>
#endif

static void
bsplit(acb_t A, acb_t B, acb_t C, acb_t D,
    const acb_t b, const acb_t z, slong n0, slong n1, slong prec)
{
    if (n1 - n0 == 1)
    {
        acb_zero(A);
        acb_one(B);
        acb_neg(C, b);
        acb_add_si(C, C, 2 - n0, prec);
        acb_mul_si(C, C, n0 - 1, prec);
        acb_sub(D, z, b, prec);
        acb_add_si(D, D, 2 - 2 * n0, prec);
    }
    else
    {
        slong m;
        acb_t T, A2, B2, C2, D2;

        acb_init(T);
        acb_init(A2);
        acb_init(B2);
        acb_init(C2);
        acb_init(D2);

        m = n0 + (n1 - n0) / 2;

        bsplit(A, B, C, D, b, z, n0, m, prec);
        bsplit(A2, B2, C2, D2, b, z, m, n1, prec);

        acb_set(T, A);
        acb_mul(A, A, A2, prec);
        acb_addmul(A, B2, C, prec);
        acb_mul(C, C, D2, prec);
        acb_addmul(C, C2, T, prec);

        acb_set(T, B);
        acb_mul(B, B, A2, prec);
        acb_addmul(B, B2, D, prec);
        acb_mul(D, D, D2, prec);
        acb_addmul(D, C2, T, prec);

        acb_clear(T);
        acb_clear(A2);
        acb_clear(B2);
        acb_clear(C2);
        acb_clear(D2);
    }
}

void
acb_hypgeom_u_si_rec(acb_t res, slong a, const acb_t b, const acb_t z, slong prec)
{
    slong k;
    acb_t u0, u1, t;

    if (a > 0)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    if (a == 0)
    {
        acb_one(res);
        return;
    }
    else if (a == -1)
    {
        acb_sub(res, z, b, prec);
        return;
    }

    /* special-case U(-n, -n+1, z) = z^n */
    if (acb_equal_si(b, a + 1))
    {
        acb_pow_si(res, z, -a, prec);
        return;
    }

    acb_init(u0);
    acb_init(u1);
    acb_init(t);

    acb_one(u0);
    acb_sub(u1, z, b, prec);

    if (-a < 20)
    {
        for (k = 2; k <= -a; k++)
        {
            acb_neg(t, b);
            acb_add_si(t, t, 2 - k, prec);
            acb_mul_si(t, t, k - 1, prec);
            acb_mul(u0, u0, t, prec);

            acb_sub(t, z, b, prec);
            acb_add_si(t, t, 2 - 2 * k, prec);
            acb_addmul(u0, u1, t, prec);

            acb_swap(u0, u1);
        }

        acb_set(res, u1);
    }
    else
    {
        acb_t A, B, C, D;

        acb_init(A);
        acb_init(B);
        acb_init(C);
        acb_init(D);

        bsplit(A, B, C, D, b, z, 2, -a + 1, prec);

        acb_sub(A, z, b, prec);
        acb_mul(D, D, A, prec);
        acb_add(res, C, D, prec);

        acb_clear(A);
        acb_clear(B);
        acb_clear(C);
        acb_clear(D);
    }

    acb_clear(u0);
    acb_clear(u1);
    acb_clear(t);
}

void
acb_hypgeom_u_1f1_series(acb_poly_t res,
    const acb_poly_t a, const acb_poly_t b, const acb_poly_t z,
    slong len, slong prec)
{
    acb_poly_t s, u, A, B;
    acb_poly_struct aa[3];
    arb_t c;
    slong wlen;
    int singular;

    acb_poly_init(s);
    acb_poly_init(u);
    acb_poly_init(A);
    acb_poly_init(B);
    acb_poly_init(aa + 0);
    acb_poly_init(aa + 1);
    acb_poly_init(aa + 2);
    arb_init(c);

    singular = (b->length == 0) || acb_is_int(b->coeffs);
    wlen = len + (singular != 0);

    /* A = rgamma(a-b+1) * 1F~1(a,b,z) */
    acb_poly_sub(u, a, b, prec);
    acb_poly_add_si(u, u, 1, prec);
    acb_poly_rgamma_series(A, u, wlen, prec);

    /* todo: handle a = 1 efficiently */
    acb_poly_set(aa, a);
    acb_poly_set(aa + 1, b);
    acb_poly_one(aa + 2);
    acb_hypgeom_pfq_series_direct(s, aa, 1, aa + 1, 2, z, 1, -1, wlen, prec);
    acb_poly_mullow(A, A, s, wlen, prec);

    /* B = rgamma(a) * 1F~1(a-b+1,2-b,z) * z^(1-b) */
    acb_poly_set(aa, u);
    acb_poly_add_si(aa + 1, b, -2, prec);
    acb_poly_neg(aa + 1, aa + 1);
    acb_hypgeom_pfq_series_direct(s, aa, 1, aa + 1, 2, z, 1, -1, wlen, prec);
    acb_poly_rgamma_series(B, a, wlen, prec);
    acb_poly_mullow(B, B, s, wlen, prec);

    acb_poly_add_si(u, b, -1, prec);
    acb_poly_neg(u, u);
    acb_poly_pow_series(s, z, u, wlen, prec);
    acb_poly_mullow(B, B, s, wlen, prec);

    acb_poly_sub(A, A, B, prec);

    /* multiply by pi csc(pi b) */
    acb_poly_sin_pi_series(B, b, wlen, prec);

    if (singular)
    {
        acb_poly_shift_right(A, A, 1);
        acb_poly_shift_right(B, B, 1);
    }

    acb_poly_div_series(res, A, B, len, prec);

    arb_const_pi(c, prec);
    _acb_vec_scalar_mul_arb(res->coeffs, res->coeffs, res->length, c, prec);

    acb_poly_clear(s);
    acb_poly_clear(u);
    acb_poly_clear(A);
    acb_poly_clear(B);
    acb_poly_clear(aa + 0);
    acb_poly_clear(aa + 1);
    acb_poly_clear(aa + 2);
    arb_clear(c);
}

void
acb_hypgeom_u_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec)
{
    if (acb_is_int(b))
    {
        acb_poly_t aa, bb, zz;

        acb_poly_init(aa);
        acb_poly_init(bb);
        acb_poly_init(zz);

        acb_poly_set_acb(aa, a);
        acb_poly_set_coeff_acb(bb, 0, b);
        acb_poly_set_coeff_si(bb, 1, 1);
        acb_poly_set_acb(zz, z);

        acb_hypgeom_u_1f1_series(zz, aa, bb, zz, 1, prec);
        acb_poly_get_coeff_acb(res, zz, 0);

        acb_poly_clear(aa);
        acb_poly_clear(bb);
        acb_poly_clear(zz);
    }
    else
    {
        acb_t t, u, v;
        acb_struct aa[3];

        acb_init(t);
        acb_init(u);
        acb_init(v);
        acb_init(aa + 0);
        acb_init(aa + 1);
        acb_init(aa + 2);

        acb_set(aa, a);
        acb_set(aa + 1, b);
        acb_one(aa + 2);
        acb_hypgeom_pfq_direct(u, aa, 1, aa + 1, 2, z, -1, prec);
        acb_sub(aa, a, b, prec);
        acb_add_ui(aa, aa, 1, prec);
        acb_sub_ui(aa + 1, b, 2, prec);
        acb_neg(aa + 1, aa + 1);
        acb_hypgeom_pfq_direct(v, aa, 1, aa + 1, 2, z, -1, prec);

        acb_sub_ui(aa + 1, b, 1, prec);

        /* rgamma(a-b+1) * gamma(1-b) * u */
        acb_rgamma(t, aa, prec);
        acb_mul(u, u, t, prec);
        acb_neg(t, aa + 1);
        acb_gamma(t, t, prec);
        acb_mul(u, u, t, prec);

        /* rgamma(a) * gamma(b-1) * z^(1-b) * v */
        acb_rgamma(t, a, prec);
        acb_mul(v, v, t, prec);
        acb_gamma(t, aa + 1, prec);
        acb_mul(v, v, t, prec);
        acb_neg(t, aa + 1);
        acb_pow(t, z, t, prec);
        acb_mul(v, v, t, prec);

        acb_add(res, u, v, prec);

        acb_clear(t);
        acb_clear(u);
        acb_clear(v);
        acb_clear(aa + 0);
        acb_clear(aa + 1);
        acb_clear(aa + 2);
    }
}

void
acb_hypgeom_u_choose(int * asymp, slong * wp,
    const acb_t a, const acb_t b, const acb_t z, slong prec)
{
    double x, y, t, cancellation;
    double input_accuracy, direct_accuracy, asymp_accuracy;

    *asymp = 0;
    *wp = prec;

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

    asymp_accuracy = sqrt(x * x + y * y) * 1.44269504088896;

    /* The Kummer transformation gives less cancellation with the 1F1 series.
    if (x < 0.0)
    {
        *kummer = 1;
        x = -x;
    } */

    if (asymp_accuracy >= prec)
    {
        *asymp = 1;
        *wp = FLINT_MAX(2, FLINT_MIN(input_accuracy + 20, prec));
        return;
    }

    /* Assume U ~ 1, M ~ exp(|z|)  (there is cancellation both in the
       evaluation of M and in the linear combination) -- a better estimate
       would account for a, b. */
    cancellation = sqrt(x * x + y * y) * 1.44269504088896 + 5;

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
acb_hypgeom_u_nointegration(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec)
{
    acb_t t;
    arf_srcptr av, tv;

    av = arb_midref(acb_realref(a));

    /* Handle small polynomial cases without divisions. */
    /* todo: should incorporate a -> 1+a-b transformation, also... */
    if (acb_is_int(a) && arf_sgn(av) <= 0)
    {
        if (arf_cmpabs_ui(av, 30) < 0 ||
            (arf_cmpabs_ui(av, prec) < 0 &&
                ((acb_bits(b) < 0.1 * prec && acb_bits(z) < 0.1 * prec)
                    || acb_contains_zero(z))))
        {
            acb_hypgeom_u_si_rec(res, arf_get_si(av, ARF_RND_DOWN), b, z, prec);
            return;
        }
    }

    acb_init(t);
    acb_sub(t, a, b, prec);
    acb_add_ui(t, t, 1, prec);
    tv = arb_midref(acb_realref(t));

    /* todo: combine these conditions with the code below */
    if ((acb_is_int(a) && arf_sgn(av) <= 0) ||
        (acb_is_int(t) && arf_sgn(tv) <= 0) ||
        acb_hypgeom_u_use_asymp(z, prec))
    {
        acb_neg(t, a);
        acb_pow(t, z, t, prec);
        acb_hypgeom_u_asymp(res, a, b, z, -1, prec);
        acb_mul(res, res, t, prec);
    }
    else
    {
        slong wp;
        int asymp;

        acb_hypgeom_u_choose(&asymp, &wp, a, b, z, prec);

        if (asymp)
        {
            acb_neg(t, a);
            acb_pow(t, z, t, prec);
            acb_hypgeom_u_asymp(res, a, b, z, -1, wp);
            acb_mul(res, res, t, prec);
        }
        else
        {
            acb_hypgeom_u_1f1(res, a, b, z, wp);
            acb_set_round(res, res, prec);
        }
    }

    acb_clear(t);
}

void
acb_hypgeom_u(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec)
{
    acb_t res2;
    slong acc, max, t;

    acb_init(res2);

    acb_hypgeom_u_nointegration(res2, a, b, z, prec);

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
                arb_is_positive(acb_realref(z)) &&
                arf_cmpabs_2exp_si(arb_midref(acb_realref(a)), 60) < 0 &&
                arf_cmpabs_2exp_si(arb_midref(acb_realref(b)), 60) < 0 &&
                arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 60) < 0)
            {
                arb_hypgeom_u_integration(acb_realref(res),
                        acb_realref(a), acb_realref(b), acb_realref(z), prec);
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
