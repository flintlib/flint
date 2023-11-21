/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz_vec.h"
#include "acb.h"
#include "arb_hypgeom.h"

#define UNROLL 4

static void
asymp_series(acb_t res, ulong n, acb_srcptr xpow, slong m, slong K, slong prec)
{
    slong j, k, khi, klo, u, r;
    fmpz * c;
    acb_t s;

    acb_init(s);
    c = _fmpz_vec_init(UNROLL + 1);

    k = K - 1;
    while (k >= 1)
    {
        u = FLINT_MIN(UNROLL, k);

        khi = k;
        klo = k - u + 1;

        for (j = klo; j <= khi; j++)
        {
            ulong aa = (2 * j - 1);
            ulong bb = (2 * j - 1);

            if (j == klo)
                fmpz_ui_mul_ui(c + khi - j, aa, bb);
            else
                fmpz_mul2_uiui(c + khi - j, c + khi - j + 1, aa, bb);
        }

        for (j = khi; j >= klo; j--)
        {
            ulong aa = (j);
            ulong bb = (2 * j + 2 * n + 1);

            if (n < UWORD_MAX / 8)
            {
                if (j == khi)
                {
                    fmpz_ui_mul_ui(c + u, aa, bb);
                }
                else
                {
                    fmpz_mul(c + khi - j, c + khi - j, c + u);
                    fmpz_mul2_uiui(c + u, c + u, aa, bb);
                }
            }
            else
            {
                fmpz_t t;
                fmpz_init(t);

                fmpz_set_ui(t, n);
                fmpz_add_ui(t, t, j);
                fmpz_mul_2exp(t, t, 1);
                fmpz_add_ui(t, t, 1);

                if (j == khi)
                {
                    fmpz_mul_ui(c + u, t, aa);
                }
                else
                {
                    fmpz_mul(c + khi - j, c + khi - j, c + u);
                    fmpz_mul_ui(t, t, aa);
                    fmpz_mul(c + u, c + u, t);
                }

                fmpz_clear(t);
            }
        }

        while (k >= klo)
        {
            r = k % m;

            if (k == khi)
            {
                acb_add(s, s, xpow + r, prec);
                acb_mul_fmpz(s, s, c + khi - k, prec);
            }
            else if (r == 0)
            {
                acb_add_fmpz(s, s, c + khi - k, prec);
            }
            else
            {
                acb_addmul_fmpz(s, xpow + r, c + khi - k, prec);
            }

            if (r == 0 && k != 0)
                acb_mul(s, s, xpow + m, prec);

            k--;
        }

        acb_div_fmpz(s, s, c + u, prec);
    }

    acb_add_ui(res, s, 1, prec);

    acb_clear(s);
    _fmpz_vec_clear(c, UNROLL + 1);
}


/* error: 0.50795 / (y^K sqrt(ny)) * [K! n! / (2^K (n+K)!)] */
/*                                   [K! / (2n)^K]          */
static void
_arb_hypgeom_legendre_p_ui_asymp_error(mag_t res, ulong n, const mag_t y, slong K)
{
    mag_t t, u;
    mag_init(t);
    mag_init(u);

    /* t = K! / (y^K sqrt(ny)) */
    mag_mul_ui_lower(t, y, n);
    mag_sqrt_lower(t, t);
    mag_pow_ui_lower(u, y, K);
    mag_mul_lower(t, t, u);
    mag_fac_ui(u, K);
    mag_div(t, u, t);

    if (K < n / 16)
    {
        /* (2n)^K */
        mag_set_ui_lower(u, n);
        mag_mul_2exp_si(u, u, 1);
        mag_pow_ui_lower(u, u, K);
        mag_div(t, t, u);
    }
    else
    {
        /* n! */
        mag_fac_ui(u, n);
        mag_mul(t, t, u);
        /* (n+K)! */
        mag_rfac_ui(u, n + K);
        mag_mul(t, t, u);
        /* 2^K */
        mag_mul_2exp_si(t, t, -K);
    }

    mag_mul_ui(t, t, 131);
    mag_mul_2exp_si(t, t, -8);

    mag_set(res, t);

    mag_clear(t);
    mag_clear(u);
}

int
arb_abs_le_ui(const arb_t x, ulong n)
{
    arf_struct u[3];
    arf_t t;
    int res;

    if (!arb_is_finite(x) || arf_cmpabs_ui(arb_midref(x), n) > 0)
        return 0;

    if (arb_is_exact(x))
        return 1;

    if (arf_sgn(arb_midref(x)) >= 0)
        arf_init_set_shallow(u + 0, arb_midref(x));
    else
        arf_init_neg_shallow(u + 0, arb_midref(x));

    arf_init_set_mag_shallow(u + 1, arb_radref(x));

    arf_init(u + 2);   /* no need to free */
    arf_set_ui(u + 2, n);
    arf_neg(u + 2, u + 2);

    arf_init(t);
    arf_sum(t, u, 3, MAG_BITS, ARF_RND_DOWN);
    res = (arf_sgn(t) < 0);
    arf_clear(t);

    return res;
}

void
_arb_hypgeom_legendre_p_ui_asymp(arb_t res, ulong n, const arb_t x,
    const arb_t y, acb_srcptr w4pow, const arb_t binom, slong m, slong K, slong prec)
{
    arb_t t, u;
    acb_t s, z;
    fmpz_t e;
    mag_t err;

    arb_init(t);
    arb_init(u);
    acb_init(s);
    acb_init(z);
    mag_init(err);
    fmpz_init(e);

    /* u = n + 1/2 */
    arb_set_d(u, 0.5);
    arb_add_ui(u, u, n, prec);

    arb_get_mag_lower(err, y);
    _arb_hypgeom_legendre_p_ui_asymp_error(err, n, err, K);

    /* z = (x + yi)^(n+0.5) * (1-i) */
    if (n < 256 || prec > 2000)
    {
        arb_set(acb_realref(z), x);
        arb_set(acb_imagref(z), y);
        acb_pow_arb(z, z, u, prec);
    }
    else
    {
        arb_atan2(t, y, x, prec);
        arb_mul(t, t, u, prec);
        arb_sin_cos(acb_imagref(z), acb_realref(z), t, prec);
    }
    arb_one(acb_realref(s));
    arb_set_si(acb_imagref(s), -1);
    acb_mul(z, z, s, prec);

    /* main series */
    asymp_series(s, n, w4pow, m, K, prec);

    /* we will use Re(z * s) */
    acb_mul(z, z, s, prec);

    /* prefactor: t = 4^n / (pi * sqrt(y) * (n+0.5) binomial(2n,n)) */
    arb_mul(t, binom, u, prec);
    arb_sqrt(u, y, prec);
    arb_mul(t, t, u, prec);
    arb_const_pi(u, prec);
    arb_mul(t, t, u, prec);

    arb_div(res, acb_realref(z), t, prec);

    fmpz_set_ui(e, n);
    arb_mul_2exp_fmpz(res, res, e);
    arb_mul_2exp_fmpz(res, res, e);

    arb_add_error_mag(res, err);

    arb_clear(t);
    arb_clear(u);
    acb_clear(s);
    acb_clear(z);
    mag_clear(err);
    fmpz_clear(e);
}

void
arb_hypgeom_legendre_p_ui_asymp(arb_t res, arb_t res2, ulong n, const arb_t x, slong K, slong prec)
{
    arb_t y, binom;
    acb_t w;
    acb_ptr w4pow;
    slong m;

    if (n == 0)
    {
        if (res != NULL) arb_one(res);
        if (res2 != NULL) arb_zero(res2);
        return;
    }

    if (!arb_abs_le_ui(x, 1))
    {
        if (res != NULL) arb_indeterminate(res);
        if (res2 != NULL) arb_indeterminate(res2);
        return;
    }

    K = FLINT_MAX(K, 1);

    if (res2 != NULL)
        m = n_sqrt(2 * K);
    else
        m = n_sqrt(K);

    arb_init(y);
    arb_init(binom);
    acb_init(w);
    w4pow = _acb_vec_init(m + 1);

    /* y = sqrt(1-x^2) */
    arb_one(y);
    arb_submul(y, x, x, 2 * prec);
    arb_sqrt(y, y, prec);

    /* w = 1 - (x/y)i */
    arb_one(acb_realref(w));
    arb_div(acb_imagref(w), x, y, prec);
    arb_neg(acb_imagref(w), acb_imagref(w));

    acb_mul_2exp_si(w, w, -2);
    _acb_vec_set_powers(w4pow, w, m + 1, prec);

    /* binomial(2n,n) */
    arb_hypgeom_central_bin_ui(binom, n, prec);

    if (res2 == NULL)
    {
        _arb_hypgeom_legendre_p_ui_asymp(res, n, x, y, w4pow, binom, m, K, prec);
    }
    else
    {
        arb_t t, u, v;

        arb_init(t);
        arb_init(u);
        arb_init(v);

        _arb_hypgeom_legendre_p_ui_asymp(t, n, x, y, w4pow, binom, m, K, prec);

        /* recurrence for central binomial */
        arb_mul_ui(binom, binom, n, prec);
        arb_set_ui(u, n);
        arb_mul_2exp_si(u, u, 2);
        arb_sub_ui(u, u, 2, prec);
        arb_div(binom, binom, u, prec);

        _arb_hypgeom_legendre_p_ui_asymp(u, n - 1, x, y, w4pow, binom, m, K, prec);

        /* P' = n (P[n-1] - x P) / (1 - x^2) */
        arb_submul(u, t, x, prec);
        arb_mul(v, x, x, 2 * prec);
        arb_neg(v, v);
        arb_add_ui(v, v, 1, prec);
        arb_div(u, u, v, prec);
        arb_mul_ui(res2, u, n, prec);
        if (res != NULL)
            arb_set(res, t);

        arb_clear(t);
        arb_clear(u);
        arb_clear(v);
    }

    arb_clear(y);
    arb_clear(binom);
    acb_clear(w);
    _acb_vec_clear(w4pow, m + 1);
}

