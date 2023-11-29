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

int
acb_hypgeom_u_asymp_determine_region(const mag_t r,
    const mag_t zlo, const acb_t z);

void
acb_hypgeom_gamma_upper_asymp(acb_t res, const acb_t s,
        const acb_t z, int regularized, slong prec)
{
    acb_t t, u;

    acb_init(t);
    acb_init(u);

    acb_sub_ui(t, s, 1, prec);
    acb_neg(t, t);

    acb_hypgeom_u_asymp(u, t, t, z, -1, prec);

    if (regularized == 2)
    {
        acb_div(u, u, z, prec);
    }
    else
    {
        acb_neg(t, t);
        acb_pow(t, z, t, prec);
        acb_mul(u, u, t, prec);

        if (regularized == 1)
        {
            acb_rgamma(t, s, prec);
            acb_mul(u, u, t, prec);
        }
    }

    acb_neg(t, z);
    acb_exp(t, t, prec);
    acb_mul(res, t, u, prec);

    acb_clear(t);
    acb_clear(u);
}

void
acb_hypgeom_gamma_upper_1f1a(acb_t res, const acb_t s,
        const acb_t z, int regularized, slong prec)
{
    acb_t a, t, w;
    acb_struct b[2];

    acb_init(a);
    acb_init(b);
    acb_init(b + 1);
    acb_init(t);
    acb_init(w);

    acb_set(a, s);
    acb_add_ui(b, s, 1, prec);
    acb_one(b + 1);
    acb_neg(w, z);

    /* t = 1F1(s, s+1, -z) / s */
    acb_hypgeom_pfq_direct(t, a, 1, b, 2, w, -1, prec);
    acb_div(t, t, s, prec);

    if (regularized == 2)
    {
        acb_neg(a, s);
        acb_pow(a, z, a, prec);
        acb_gamma(b, s, prec);
        acb_mul(b, b, a, prec);
        acb_sub(res, b, t, prec);
    }
    else
    {
        acb_pow(a, z, s, prec);
        acb_mul(t, t, a, prec);

        if (regularized == 1)
        {
            acb_rgamma(a, s, prec);
            acb_mul(t, t, a, prec);
            acb_sub_ui(res, t, 1, prec);
            acb_neg(res, res);
        }
        else
        {
            acb_gamma(a, s, prec);
            acb_sub(res, a, t, prec);
        }
    }

    acb_clear(a);
    acb_clear(b);
    acb_clear(b + 1);
    acb_clear(t);
    acb_clear(w);
}

void
acb_hypgeom_gamma_upper_1f1b(acb_t res, const acb_t s,
        const acb_t z, int regularized, slong prec)
{
    acb_t a, b, t;

    acb_init(a);
    acb_init(b);
    acb_init(t);

    acb_add_ui(b, s, 1, prec);
    acb_hypgeom_pfq_direct(t, NULL, 0, b, 1, z, -1, prec);
    acb_div(t, t, s, prec);

    acb_neg(a, z);
    acb_exp(a, a, prec);
    acb_mul(t, t, a, prec);

    if (regularized == 2)
    {
        acb_neg(a, s);
        acb_pow(a, z, a, prec);
        acb_gamma(b, s, prec);
        acb_mul(b, b, a, prec);
        acb_sub(res, b, t, prec);
    }
    else
    {
        acb_pow(a, z, s, prec);
        acb_mul(t, t, a, prec);

        if (regularized == 1)
        {
            acb_rgamma(a, s, prec);
            acb_mul(t, t, a, prec);
            acb_sub_ui(res, t, 1, prec);
            acb_neg(res, res);
        }
        else
        {
            acb_gamma(a, s, prec);
            acb_sub(res, a, t, prec);
        }
    }

    acb_clear(a);
    acb_clear(b);
    acb_clear(t);
}

/* requires n <= 0 */
void
acb_hypgeom_gamma_upper_singular(acb_t res, slong s, const acb_t z, int regularized, slong prec)
{
    acb_t A, B, C, t, u;
    acb_struct a[2];
    acb_struct b[2];
    arb_t f;
    slong n;

    if (regularized == 1)
    {
        acb_zero(res);
        return;
    }

    n = -s;

    acb_init(A);
    acb_init(B);
    acb_init(C);
    acb_init(t);
    acb_init(u);
    arb_init(f);
    acb_init(a + 0);
    acb_init(a + 1);
    acb_init(b + 0);
    acb_init(b + 1);

    arb_fac_ui(f, n, prec);

    /* (-1)^n (psi(n+1) - log(z)) / n! */
    acb_set_ui(A, n + 1);
    acb_digamma(A, A, prec);

    acb_log(t, z, prec);
    acb_sub(A, A, t, prec);
    acb_div_arb(A, A, f, prec);
    if (n % 2)
        acb_neg(A, A);

    /* (-1)^n z 2F2(1,1;2,2+n;-z) / (n+1)! */
    acb_set_si(a, 1);
    acb_set_si(b, 2);
    acb_set_si(b + 1, 2 + n);
    acb_neg(t, z);
    acb_hypgeom_pfq_direct(B, a, 1, b, 2, t, -1, prec);
    acb_mul(B, B, z, prec);
    arb_mul_ui(f, f, n + 1, prec);
    acb_div_arb(B, B, f, prec);
    if (n % 2)
        acb_neg(B, B);

    /* -sum((-z)^k / ((k-n) k!), k=0...n-1) */
    acb_set_si(a, -n);
    acb_set_si(b, 1 - n);
    acb_set_si(b + 1, 1);
    acb_neg(t, z);

    if (n == 0)
    {
        acb_zero(C);
    }
    else
    {
        acb_hypgeom_pfq_sum(C, u, a, 1, b, 2, t, n, prec);
        acb_div_ui(C, C, n, prec);
    }

    if (regularized == 2)
    {
        acb_pow_si(t, z, n, prec);
        acb_mul(A, A, t, prec);
        acb_mul(B, B, t, prec);
    }
    else
    {
        acb_pow_si(t, z, -n, prec);
        acb_mul(C, C, t, prec);
    }

    acb_add(res, A, B, prec);
    acb_add(res, res, C, prec);

    acb_clear(A);
    acb_clear(B);
    acb_clear(C);
    acb_clear(t);
    acb_clear(u);
    arb_clear(f);
    acb_clear(a + 0);
    acb_clear(a + 1);
    acb_clear(b + 0);
    acb_clear(b + 1);
}

static int
_determine_region(const acb_t s, const acb_t z)
{
    int R;
    mag_t r, zlo;
    acb_t t;

    mag_init(r);
    mag_init(zlo);
    acb_init(t);

    /* lower bound for |z| */
    acb_get_mag_lower(zlo, z);

    /* upper bound for r = |s - 1| */
    acb_sub_ui(t, s, 1, MAG_BITS);
    acb_get_mag(r, t);

    R = acb_hypgeom_u_asymp_determine_region(r, zlo, z);

    mag_clear(r);
    mag_clear(zlo);
    acb_clear(t);

    return R;
}

/* Returns 1 if it can be determined that a^n > b^n + c^n.
 * Requires 1 <= n <= WORD_MAX.
 * If n == WORD_MAX, the infinity norm a > max(b, c) is compared instead. */
int
_mag_gt_norm_ui(const mag_t a, const mag_t b, const mag_t c, ulong n)
{
    int result;
    result = 0;
    if (n < 1)
    {
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);
    }
    else if (mag_is_zero(a))
    {
        result = 0;
    }
    else if (mag_is_zero(b))
    {
        result = mag_cmp(a, c) > 0;
    }
    else if (mag_is_zero(c))
    {
        result = mag_cmp(a, b) > 0;
    }
    else if (n == WORD_MAX)
    {
        result = mag_cmp(a, b) > 0 && mag_cmp(a, c) > 0;
    }
    else if (n == 1)
    {
        mag_t sum;
        mag_init(sum);
        mag_add(sum, b, c);
        result = mag_cmp(a, sum) > 0;
        mag_clear(sum);
    }
    else if (_mag_gt_norm_ui(a, b, c, 1))
    {
        result = 1;
    }
    else if (!_mag_gt_norm_ui(a, b, c, WORD_MAX))
    {
        result = 0;
    }
    else
    {
        mag_t u, v, w, sum;
        mag_init(u);
        mag_init(v);
        mag_init(w);
        mag_init(sum);
        mag_pow_ui_lower(u, a, n);
        mag_pow_ui(v, b, n);
        mag_pow_ui(w, c, n);
        mag_add(sum, v, w);
        result = mag_cmp(u, sum) > 0;
        mag_clear(u);
        mag_clear(v);
        mag_clear(w);
        mag_clear(sum);
    }
    return result;
}

/* Given nonnegative real s, x, and prec,
 * returns 1 if the asymptotic expansion of the
 * hypergeometric U function should be used for upper incomplete gamma.
 * It returns 1 if x > 4 and (x-c)^8 > (a*s)^8 + (b*prec)^8.
 * Careful mag rounding is not used because this is just a heuristic. */
static int
_nonnegative_real_use_asymp(const mag_t s, const mag_t x, slong prec)
{
    int result;
    result = 0;
    if (mag_cmp_2exp_si(x, 2) > 0)
    {
        mag_t a, b, c, u, v, w;
        mag_init(a);
        mag_init(b);
        mag_init(c);
        mag_init(u);
        mag_init(v);
        mag_init(w);
        mag_set_d(a, 1.029287542);
        mag_set_d(b, 0.3319411658);
        mag_set_d(c, 2.391097143);
        mag_sub(u, x, c);
        mag_mul(v, a, s);
        mag_mul_ui(w, b, FLINT_MAX(0, prec));
        result = _mag_gt_norm_ui(u, v, w, 8);
        mag_clear(a);
        mag_clear(b);
        mag_clear(c);
        mag_clear(u);
        mag_clear(v);
        mag_clear(w);
    }
    return result;
}

static int
_acb_is_nonnegative_real(const acb_t z)
{
    return arb_is_zero(acb_imagref(z)) && arb_is_nonnegative(acb_realref(z));
}

void
acb_hypgeom_gamma_upper_nointegration(acb_t res, const acb_t s, const acb_t z, int regularized, slong prec)
{
    if (!acb_is_finite(s) || !acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (acb_is_zero(z))
    {
        if (regularized == 2)
        {
            if (arb_is_negative(acb_realref(s)))
            {
                acb_inv(res, s, prec);
                acb_neg(res, res);
            }
            else
            {
                acb_indeterminate(res);
            }
        }
        else
        {
            if (arb_is_positive(acb_realref(s)))
            {
                if (regularized == 1)
                    acb_one(res);
                else
                    acb_gamma(res, s, prec);
            }
            else
            {
                acb_indeterminate(res);
            }
        }
    }
    else
    {
        slong n = WORD_MAX;

        if (acb_is_int(s))
        {
            if (regularized == 1 && arf_sgn(arb_midref(acb_realref(s))) <= 0)
            {
                acb_zero(res);
                return;
            }

            if (arf_cmpabs_2exp_si(arb_midref(acb_realref(s)), 30) < 0)
            {
                n = arf_get_si(arb_midref(acb_realref(s)), ARF_RND_DOWN);
            }
        }

        if (n >= 1 && n <= 3)
        {
            acb_t t, u;

            acb_init(t);
            acb_init(u);

            if (regularized == 2)
                acb_pow_si(u, z, n, prec);

            if (n == 1)
            {
                acb_neg(res, z);
                acb_exp(res, res, prec);
            }
            else if (n == 2)
            {
                acb_add_ui(t, z, 1, prec);
                acb_neg(res, z);
                acb_exp(res, res, prec);
                acb_mul(res, res, t, prec);
            }
            else if (n == 3)
            {
                acb_add_ui(t, z, 2, prec);
                acb_mul(t, t, z, prec);
                acb_add_ui(t, t, 2, prec);
                acb_neg(res, z);
                acb_exp(res, res, prec);
                acb_mul(res, res, t, prec);
            }

            if (regularized == 2)
                acb_div(res, res, u, prec);
            else if (regularized == 1 && n == 3)
                acb_mul_2exp_si(res, res, -1);

            acb_clear(t);
            acb_clear(u);

            return;
        }

        if (_acb_is_nonnegative_real(s) && _acb_is_nonnegative_real(z))
        {
            int result;
            mag_t ms, mx;
            mag_init(ms);
            mag_init(mx);
            arb_get_mag(ms, acb_realref(s));
            arb_get_mag(mx, acb_realref(z));
            result = _nonnegative_real_use_asymp(ms, mx, prec);
            mag_clear(ms);
            mag_clear(mx);
            if (result)
            {
                acb_hypgeom_gamma_upper_asymp(res, s, z, regularized, prec);
                return;
            }
        }
        else if (acb_hypgeom_u_use_asymp(z, prec) &&
                 ((0 < n && n < WORD_MAX) || _determine_region(s, z)))
        {
            acb_hypgeom_gamma_upper_asymp(res, s, z, regularized, prec);
            return;
        }

        if (n <= 0 && n > -10 * prec)
        {
            acb_hypgeom_gamma_upper_singular(res, n, z, regularized, prec);
            return;
        }

        if (arf_sgn(arb_midref(acb_realref(z))) > 0)
            acb_hypgeom_gamma_upper_1f1b(res, s, z, regularized, prec);
        else
            acb_hypgeom_gamma_upper_1f1a(res, s, z, regularized, prec);
    }
}

void
acb_hypgeom_gamma_upper(acb_t res, const acb_t s, const acb_t z, int regularized, slong prec)
{
    acb_t res2;
    slong acc, max, t;

    acb_init(res2);

    acb_hypgeom_gamma_upper_nointegration(res2, s, z, regularized, prec);

    acc = acb_rel_accuracy_bits(res2);

    if (acc < 0.5 * prec)
    {
        max = prec;
        t = acb_rel_accuracy_bits(z);
        max = FLINT_MIN(max, t);
        t = acb_rel_accuracy_bits(s);
        max = FLINT_MIN(max, t);

        if (max > 2 && acc < 0.5 * max)
        {
            if (acb_is_real(s) && acb_is_real(z) &&
                arf_cmp_2exp_si(arb_midref(acb_realref(z)), -16) > 0 &&
                arf_cmpabs_2exp_si(arb_midref(acb_realref(s)), 60) < 0 &&
                arf_cmpabs_2exp_si(arb_midref(acb_realref(z)), 60) < 0)
            {
                arb_hypgeom_gamma_upper_integration(acb_realref(res),
                        acb_realref(s), acb_realref(z), regularized, prec);
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
