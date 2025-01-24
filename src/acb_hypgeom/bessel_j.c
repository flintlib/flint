/*
    Copyright (C) 2014-2015, 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "acb.h"
#include "acb_hypgeom.h"

/*
Bound |J^{(d)}_nu(z)|, or |J_nu(z)| itself with d = 0.
*/
void
acb_hypgeom_bessel_j_deriv_bound(mag_t res, const acb_t nu, const acb_t z, ulong d)
{
    if (!acb_is_int(nu))
    {
        /* Not yet implemented. */
        /* Things we could use:

           - For complex z, complex n, integer d >= 0:
                |J^{(d)}_n(z)| <= max(|J_{n-d}(z)|, ..., |J_{n+d}(z)|).
            - Asymptotic expansion (general case)
                Bounds are messy and not useful when |nu| is large.
        */
        mag_inf(res);
        return;
    }

    mag_t bound, t, u, x, y, zlow;

    mag_init(bound);
    mag_init(t);
    mag_init(u);
    mag_init(x);
    mag_init(y);
    mag_init(zlow);

    arb_get_mag(x, acb_realref(z));
    arb_get_mag(y, acb_imagref(z));
    acb_get_mag_lower(zlow, z);

    /* |J^{(d)}_n(x)| <= exp(|im(z)|). */
    mag_exp(bound, y);

    /* |n| is something reasonable, do something less generic */
    if (arf_cmpabs_2exp_si(arb_midref(acb_realref(nu)), 30) < 0)
    {
        slong n;

        n = arf_get_si(arb_midref(acb_realref(nu)), ARF_RND_DOWN);
        n = FLINT_ABS(n);

        /*
        Small |z| bound.
        [DLMF 10.14.4] For real n >= -0.5, complex z:

            |J_n(z)| <= (|z/2|^n / n!) exp(|im(z)|)

        Analogously for d = 1:

            |J'_0(z)| <= (|z|/2) exp(|im(z)|)
            |J'_1(z)| <= (1/2) exp(|im(z)|)
            |J'_n(z)| <= (|z|^(n-1) / (2^n * (n-1)!)) exp(|im(z)|)

        TODO: generalize to higher derivatives.
        */
        if (d == 0)
        {
            mag_hypot(t, x, y);
            mag_mul_2exp_si(t, t, -1);
            mag_pow_ui(t, t, n);
            mag_rfac_ui(u, n);
            mag_mul(t, t, u);
            if (mag_cmp_2exp_si(t, 0) < 0)
                mag_mul(bound, bound, t);
        }
        else if (d == 1)
        {
            if (n == 1)
            {
                mag_mul_2exp_si(bound, bound, -1);
            }
            else
            {
                mag_hypot(t, x, y);

                if (n == 0)
                    mag_mul_2exp_si(t, t, -1);
                else
                {
                    mag_pow_ui(t, t, n - 1);
                    mag_mul_2exp_si(t, t, -n);
                    mag_rfac_ui(u, n - 1);
                    mag_mul(t, t, u);
                }

                if (mag_cmp_2exp_si(t, 0) < 0)
                    mag_mul(bound, bound, t);
            }
        }

        /*
        For integer n,

            |J^{(d)}_n(z)| <= C_{d,n} * cosh(im(z)) / sqrt(2/(pi |z|))

        for some C_{d,n} >= 1.

        Proving tight explicit bounds for C_{d,n} is an open problem:
        https://mathoverflow.net/questions/485404

        However, we can brute force some reasonable bounds for small n, d
        which are most commonly needed. TODO: extend this table.
        */
        if (n <= 16 && d <= 16)
        {
            /* sqrt(2/pi) + 0.0001 */
            double Csq2pi = 0.798;

            static const float Cbound[] = {
                1.077,   /* n = 0 */
                1.035,   /* n = 1 */
                1.089,   /* n = 2 */
                1.131,   /* n = 3 */
                1.167,   /* n = 4 */
                1.197,   /* n = 5 */
                1.223,   /* n = 6 */
                1.247,   /* n = 7 */
                1.269,   /* n = 8 */
                1.288,   /* n = 9 */
                1.306,   /* n = 10 */
                1.323,   /* n = 11 */
                1.339,   /* n = 12 */
                1.354,   /* n = 13 */
                1.368,   /* n = 14 */
                1.382,   /* n = 15 */
                1.394,   /* n = 16 */
            };

            if (d == 0)
                Csq2pi *= Cbound[n];
            else if (n == 0 && d == 1)
                Csq2pi *= Cbound[1];

            mag_rsqrt(t, zlow);
            mag_set_d(u, Csq2pi);
            mag_mul(t, t, u);
            mag_cosh(u, y);
            mag_mul(t, t, u);
            mag_min(bound, bound, t);
        }
        else if (mag_is_zero(y))
        {
            /* Landau: |J^{(d)}_n(x)| <= 0.786 |x|^(-1/3) */
            mag_inv(t, zlow);
            mag_root(t, t, 3);
            mag_mul_ui(t, t, 13186892);
            mag_mul_2exp_si(t, t, -24);
            mag_min(bound, bound, t);
        }

        /* Landau:
            |J_n(x)| <= 0.675 |n|^(-1/3)
            |J^{(d)}_n(x)| <= 0.675 max(0, |n|-d)^(-1/3)
        */
        if (mag_is_zero(y) && d < (ulong) n)
        {
            mag_set_ui_lower(t, n - d);
            mag_inv(t, t);
            mag_root(t, t, 3);
            mag_mul_ui(t, t, 11324621);
            mag_mul_2exp_si(t, t, -24);
            mag_min(bound, bound, t);
        }
    }

    mag_set(res, bound);

    mag_clear(bound);
    mag_clear(t);
    mag_clear(u);
    mag_clear(x);
    mag_clear(y);
    mag_clear(zlow);
}

static int
_acb_hypgeom_bessel_j_is_real(const acb_t nu, const acb_t z)
{
    if (acb_is_int(nu))
    {
        if (arb_is_zero(acb_imagref(z)))
            return 1;
        if (arb_is_zero(acb_realref(z)))
            return arf_is_int_2exp_si(arb_midref(acb_realref(nu)), 1);
    }

    return 0;
}

static int
_acb_hypgeom_bessel_j_is_imag(const acb_t nu, const acb_t z)
{
    if (acb_is_int(nu) && arb_is_zero(acb_realref(z)))
        return !arf_is_int_2exp_si(arb_midref(acb_realref(nu)), 1);
    else
        return 0;
}

/* Bound propagated error for |J_n(z)| when evaluating at mid(z). */
void
_acb_hypgeom_bessel_j_prop_error(mag_t re, mag_t im, const acb_t nu, const acb_t z)
{
    mag_t err, rad;

    mag_init(err);
    mag_init(rad);

    acb_hypgeom_bessel_j_deriv_bound(err, nu, z, 1);

    if (!mag_is_finite(err))
    {
        mag_inf(re);
        mag_inf(im);
    }
    else
    {
        mag_hypot(rad, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));
        mag_mul(err, err, rad);

        if (_acb_hypgeom_bessel_j_is_real(nu, z))
        {
            mag_set(re, err);
            mag_zero(im);
        }
        else if (_acb_hypgeom_bessel_j_is_imag(nu, z))
        {
            mag_zero(re);
            mag_set(im, err);
        }
        else
        {
            mag_set(re, err);
            mag_set(im, err);
        }
    }

    mag_clear(err);
    mag_clear(rad);
}

/* assumes no aliasing */
/* (+/- iz)^(-1/2-v) * z^v * exp(+/- iz) */
void
acb_hypgeom_bessel_j_asymp_prefactors_fallback(acb_t Ap, acb_t Am, acb_t C,
    const acb_t nu, const acb_t z, slong prec)
{
    acb_t t, u, v;

    acb_init(t);
    acb_init(u);
    acb_init(v);

    /* v = -1/2-nu */
    acb_one(v);
    acb_mul_2exp_si(v, v, -1);
    acb_add(v, v, nu, prec);
    acb_neg(v, v);

    acb_mul_onei(t, z);  /* t = iz */
    acb_neg(u, t);       /* u = -iz */

    /* Ap, Am = (+/- iz)^(-1/2-nu) */
    acb_pow(Ap, t, v, prec);
    acb_pow(Am, u, v, prec);

    /* Ap, Am *= exp(+/- iz) */
    acb_exp_invexp(u, v, t, prec);
    acb_mul(Ap, Ap, u, prec);
    acb_mul(Am, Am, v, prec);

    /* z^nu */
    acb_pow(t, z, nu, prec);
    acb_mul(Ap, Ap, t, prec);
    acb_mul(Am, Am, t, prec);

    /* (2 pi)^(-1/2) */
    acb_const_pi(C, prec);
    acb_mul_2exp_si(C, C, 1);
    acb_rsqrt(C, C, prec);

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
}

void
acb_hypgeom_bessel_j_asymp_prefactors(acb_t Ap, acb_t Am, acb_t C,
    const acb_t nu, const acb_t z, slong prec)
{
    if (arb_is_positive(acb_realref(z)))
    {
        acb_t t, u;

        acb_init(t);
        acb_init(u);

        /* -(2nu+1)/4 * pi + z */
        acb_mul_2exp_si(t, nu, 1);
        acb_add_ui(t, t, 1, prec);
        acb_mul_2exp_si(t, t, -2);
        acb_neg(t, t);
        acb_const_pi(u, prec);
        acb_mul(t, t, u, prec);
        acb_add(t, t, z, prec);
        acb_mul_onei(t, t);
        acb_exp_invexp(Ap, Am, t, prec);

        /* (2 pi z)^(-1/2) */
        acb_const_pi(C, prec);
        acb_mul_2exp_si(C, C, 1);
        acb_mul(C, C, z, prec);
        acb_rsqrt(C, C, prec);

        acb_clear(t);
        acb_clear(u);
        return;
    }

    acb_hypgeom_bessel_j_asymp_prefactors_fallback(Ap, Am, C, nu, z, prec);
}

void
acb_hypgeom_bessel_j_asymp(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    acb_t A1, A2, C, U1, U2, s, t, u;
    int is_real, is_imag;

    /* zero at -inf and +inf when nu is finite */
    if (acb_is_finite(nu) && !acb_is_finite(z) &&
        acb_is_real(z) && !acb_contains_zero(z))
    {
        acb_zero(res);
        return;
    }

    acb_init(A1);
    acb_init(A2);
    acb_init(C);
    acb_init(U1);
    acb_init(U2);
    acb_init(s);
    acb_init(t);
    acb_init(u);

    is_imag = 0;
    is_real = acb_is_real(nu) && acb_is_real(z)
        && (acb_is_int(nu) || arb_is_positive(acb_realref(z)));

    if (!is_real && arb_is_zero(acb_realref(z)) && acb_is_int(nu))
    {
        acb_mul_2exp_si(t, nu, -1);

        if (acb_is_int(t))
            is_real = 1;
        else
            is_imag = 1;
    }

    acb_hypgeom_bessel_j_asymp_prefactors(A1, A2, C, nu, z, prec);

    /* todo: if Ap ~ 2^a and Am = 2^b and U1 ~ U2 ~ 1, change precision? */

    if (!acb_is_finite(A1) || !acb_is_finite(A2) || !acb_is_finite(C))
    {
        acb_indeterminate(res);
    }
    else
    {
        /* s = 1/2 + nu */
        acb_one(s);
        acb_mul_2exp_si(s, s, -1);
        acb_add(s, s, nu, prec);

        /* t = 1 + 2 nu */
        acb_mul_2exp_si(t, nu, 1);
        acb_add_ui(t, t, 1, prec);

        acb_mul_onei(u, z);
        acb_mul_2exp_si(u, u, 1);
        acb_hypgeom_u_asymp(U2, s, t, u, -1, prec);
        acb_neg(u, u);
        acb_hypgeom_u_asymp(U1, s, t, u, -1, prec);

        acb_mul(res, A1, U1, prec);
        acb_addmul(res, A2, U2, prec);
        acb_mul(res, res, C, prec);

        if (is_real)
            arb_zero(acb_imagref(res));
        if (is_imag)
            arb_zero(acb_realref(res));
    }

    acb_clear(A1);
    acb_clear(A2);
    acb_clear(C);
    acb_clear(U1);
    acb_clear(U2);
    acb_clear(s);
    acb_clear(t);
    acb_clear(u);
}

static void
_acb_hypgeom_bessel_j_0f1(acb_t res, const acb_t nu, const acb_t z, slong prec, slong prec2)
{
    acb_struct b[2];
    acb_t w, c, t;

    if (acb_is_int(nu) && arb_is_negative(acb_realref(nu)))
    {
        acb_init(t);
        acb_neg(t, nu);

        _acb_hypgeom_bessel_j_0f1(res, t, z, prec, prec2);

        acb_mul_2exp_si(t, t, -1);
        if (!acb_is_int(t))
            acb_neg(res, res);

        acb_clear(t);
        return;
    }

    acb_init(b + 0);
    acb_init(b + 1);
    acb_init(w);
    acb_init(c);
    acb_init(t);

    acb_add_ui(b + 0, nu, 1, prec2);
    acb_one(b + 1);

    /* (z/2)^nu / gamma(nu+1) */
    acb_mul_2exp_si(c, z, -1);
    acb_pow(c, c, nu, prec);
    acb_rgamma(t, b + 0, prec);
    acb_mul(c, t, c, prec);

    /* -z^2/4 */
    acb_mul(w, z, z, prec2);
    acb_mul_2exp_si(w, w, -2);
    acb_neg(w, w);

    acb_hypgeom_pfq_direct(t, NULL, 0, b, 2, w, -1, prec2);

    acb_mul(res, t, c, prec);

    acb_clear(b + 0);
    acb_clear(b + 1);
    acb_clear(w);
    acb_clear(c);
    acb_clear(t);
}

void
acb_hypgeom_bessel_j_0f1(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    _acb_hypgeom_bessel_j_0f1(res, nu, z, prec, prec);
}

void
acb_hypgeom_bessel_j_transition(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    slong prec2;
    double a, b, zz;
    slong cancellation;

    /* All terms are positive, so no cancellation */
    if (arb_is_zero(acb_realref(z)) && (acb_is_int(nu) || (acb_is_real(nu) && arb_is_positive(acb_realref(nu)))))
    {
        acb_hypgeom_bessel_j_0f1(res, nu, z, prec);
        return;
    }

    a = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_DOWN);
    b = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_DOWN);
    zz = sqrt(a * a + b * b);

    /* If nu > |z|^2/4, there is no significant cancellation. */
    if (acb_is_real(nu) && arf_cmpabs_2exp_si(arb_midref(acb_realref(nu)), 20) < 64)
    {
        double nn;
        nn = arf_get_d(arb_midref(acb_realref(nu)), ARF_RND_DOWN);

        if (nn > zz * zz * 0.25)
        {
            acb_hypgeom_bessel_j_0f1(res, nu, z, prec);
            return;
        }
    }

    /* Estimate cancellation as (|x|-im(|x|)) * (1/log(2)) bits.
       TODO: estimate more accurately for large nu */
    cancellation = (zz - fabs(b)) * 1.44269504088896;
    cancellation = FLINT_MIN(cancellation, 4 * prec);
    cancellation = FLINT_MAX(cancellation, 0);

    prec2 = prec + 5 + cancellation;

    if (acb_is_exact(nu) && acb_is_exact(z))
    {
        _acb_hypgeom_bessel_j_0f1(res, nu, z, prec, prec2);
    }
    else
    {
        mag_t aerr, berr;
        acb_t t;

        mag_init(aerr);
        mag_init(berr);

        _acb_hypgeom_bessel_j_prop_error(aerr, berr, nu, z);

        /* Propagated error implemented for this nu. */
        if (mag_is_finite(aerr) && mag_is_finite(berr))
        {
            acb_init(t);
            acb_get_mid(t, z);

            _acb_hypgeom_bessel_j_0f1(res, nu, t, prec, prec2);

            arb_add_error_mag(acb_realref(res), aerr);
            arb_add_error_mag(acb_imagref(res), berr);
            acb_clear(t);
        }
        else  /* Propagated error not implemented */
        {
            slong acc = acb_rel_accuracy_bits(z);

            prec2 = FLINT_MIN(prec2, acc);
            prec2 = FLINT_MAX(prec2, prec);

            _acb_hypgeom_bessel_j_0f1(res, nu, z, prec, prec2);
        }

        mag_clear(aerr);
        mag_clear(berr);
    }
}

void
acb_hypgeom_bessel_j(acb_t res, const acb_t nu, const acb_t z, slong prec)
{
    if (!acb_is_finite(nu) || !acb_is_finite(z))
    {
        /* Some infinities */
        if (acb_is_real(nu) && acb_is_finite(nu) &&
            acb_is_real(z) && (mag_is_finite(arb_radref(acb_realref(z))) && arf_is_inf(arb_midref(acb_realref(z)))))
            acb_zero(res);
        else
            acb_indeterminate(res);
    }
    else
    {
        mag_t zmag;

        mag_init(zmag);
        acb_get_mag(zmag, z);

        if (mag_cmp_2exp_si(zmag, 3) < 0)
        {
            acb_hypgeom_bessel_j_0f1(res, nu, z, prec);
        }
        else
        {
            acb_t t;
            acb_init(t);

            /* Assuming small nu, the asymptotic series can be used roughly when
               [(1+log(2))/log(2) = 2.44269504088896] * z > p
               We are a bit more conservative and use the factor 2. */

            if (mag_cmp_2exp_si(zmag, 64) >= 0 || 2 * mag_get_d(zmag) >= prec)
                acb_hypgeom_bessel_j_asymp(t, nu, z, prec);
            else
                acb_hypgeom_bessel_j_transition(t, nu, z, prec);

#if 0
            flint_printf("T: "); acb_printd(t, 10); flint_printf("\n");
#endif

            /* If the enclosure is terrible, try bounding |J_nu(z)|
               directly. TODO: cleaner way to do this (detect
               wide intervals a priori?) */
            if (acb_rel_accuracy_bits(t) < 1)
            {
                mag_t M;
                int real, imag;

                mag_init(M);

                acb_hypgeom_bessel_j_deriv_bound(M, nu, z, 0);

#if 0
                flint_printf("M: "); mag_printd(M, 10); flint_printf("\n");
#endif

                if (mag_is_finite(M))
                {
                    real = _acb_hypgeom_bessel_j_is_real(nu, z);
                    imag = _acb_hypgeom_bessel_j_is_imag(nu, z);

                    acb_zero(res);
                    if (real)
                        arb_add_error_mag(acb_realref(res), M);
                    else if (imag)
                        arb_add_error_mag(acb_imagref(res), M);
                    else
                        acb_add_error_mag(res, M);

                    if (acb_is_finite(t))
                    {
                        arb_intersection(acb_realref(t), acb_realref(t), acb_realref(res), prec);
                        arb_intersection(acb_imagref(t), acb_imagref(t), acb_imagref(res), prec);
                    }
                    else
                    {
                        acb_set(t, res);
                    }
                }

                mag_clear(M);
            }

            acb_swap(res, t);
            acb_clear(t);
        }

        mag_clear(zmag);
    }
}
