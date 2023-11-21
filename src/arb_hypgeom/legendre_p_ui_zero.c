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
#include "arb_hypgeom.h"

#define UNROLL 4

static void
sum_rs_inner(arb_t s, arb_srcptr xpow, slong m, ulong n, slong K, slong prec)
{
    slong j, k, khi, klo, u, r;
    ulong d;
    slong sigma;
    fmpz * c;

    sigma = n % 2 ? 1 : -1;
    d = n / 2;

    arb_zero(s);
    c = _fmpz_vec_init(UNROLL + 1);

    k = K - 1;
    while (k >= 1)
    {
        u = FLINT_MIN(UNROLL, k);

        khi = k;
        klo = k - u + 1;

        for (j = klo; j <= khi; j++)
        {
            ulong aa = (d - j + 1);
            ulong bb = (2 * d + 2 * j + sigma);

            if (j == klo)
                fmpz_ui_mul_ui(c + khi - j, aa, bb);
            else
                fmpz_mul2_uiui(c + khi - j, c + khi - j + 1, aa, bb);
        }

        for (j = khi; j >= klo; j--)
        {
            ulong aa = (j);
            ulong bb = (2 * j + sigma);

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

        while (k >= klo)
        {
            r = k % m;
            if (k == khi)
            {
                arb_add(s, s, xpow + r, prec);
                arb_mul_fmpz(s, s, c + khi - k, prec);
            }
            else if (r == 0)
                arb_add_fmpz(s, s, c + khi - k, prec);
            else
                arb_addmul_fmpz(s, xpow + r, c + khi - k, prec);

            if (r == 0 && k != 0)
                arb_mul(s, s, xpow + m, prec);

            k--;
        }

        arb_div_fmpz(s, s, c + u, prec);
    }

    _fmpz_vec_clear(c, UNROLL + 1);
}

static void
_arb_hypgeom_legendre_p_ui_zero(arb_t res, ulong n,
    const arb_t x, arb_srcptr xpow, slong m, slong K, slong prec)
{
    arb_t s;
    slong d, prec2;
    mag_t u, a, t, err;

    d = n / 2;

    arb_init(s);
    mag_init(u);
    mag_init(a);
    mag_init(t);
    mag_init(err);

    K = FLINT_MIN(K, d + 1);
    sum_rs_inner(s, xpow, m, n, K, prec);

    prec2 = arb_rel_accuracy_bits(s);
    if (prec2 > prec)
        prec2 = prec;
    else
        prec2 = FLINT_MAX(0, prec2) + 20;

    arb_add_ui(s, s, 1, prec2);
    if (n % 2 == 1)
        arb_mul(s, s, x, prec2);
    arb_swap(res, s);

    if (d % 2 == 1)
        arb_neg(res, res);

    if (n % 2 == 0)
    {
        arb_hypgeom_central_bin_ui(s, d, prec2);
        arb_mul(res, res, s, prec2);
        arb_mul_2exp_si(res, res, -n);
    }
    else
    {
        arb_hypgeom_central_bin_ui(s, d + 1, prec2);
        arb_mul(res, res, s, prec2);
        arb_mul_ui(res, res, d + 1, prec2);
        arb_mul_2exp_si(res, res, -n);
    }

    if (K < d + 1)
    {
        mag_bin_uiui(err, n, d - K);
        mag_bin_uiui(t, n + 2 * K + (n % 2), n);
        mag_mul(err, err, t);
        arb_get_mag(t, x);
        mag_pow_ui(t, t, 2 * K + (n % 2));
        mag_mul(err, err, t);
        mag_mul_2exp_si(err, err, -n);

        arb_get_mag(t, x);
        mag_mul(t, t, t);
        mag_mul_ui(t, t, d - K + 1);
        mag_mul_ui(t, t, 2 * d + 2 * K + (n % 2 ? 1 : -1));
        mag_div_ui(t, t, K);
        mag_div_ui(t, t, 2 * K + (n % 2 ? 1 : -1));
        mag_geom_series(t, t, 0);
        mag_mul(err, err, t);

        arb_add_error_mag(res, err);
    }

    arb_clear(s);
    mag_clear(u);
    mag_clear(a);
    mag_clear(t);
    mag_clear(err);
}


void
arb_hypgeom_legendre_p_ui_zero(arb_t res, arb_t res2, ulong n,
                        const arb_t x, slong K, slong prec)
{
    arb_ptr xpow;
    arb_t t, u, v;
    slong m;

    if (n == 0)
    {
        if (res != NULL) arb_one(res);
        if (res2 != NULL) arb_zero(res2);
        return;
    }

    /* overflow protection */
    if (n > UWORD_MAX / 4)
    {
        if (res != NULL) arb_indeterminate(res);
        if (res2 != NULL) arb_indeterminate(res2);
    }

    if (res == NULL)
    {
        arb_init(v);
        arb_hypgeom_legendre_p_ui_zero(v, res2, n, x, K, prec);
        arb_clear(v);
        return;
    }

    arb_init(t);
    arb_init(u);
    arb_init(v);

    K = FLINT_MIN(K, n / 2 + 1);

    if (res2 != NULL)
        m = n_sqrt(2 * K);
    else
        m = n_sqrt(K);

    xpow = _arb_vec_init(m + 1);

    arb_mul(v, x, x, prec);
    arb_neg(v, v);
    _arb_vec_set_powers(xpow, v, m + 1, prec);

    /* todo: recycle prefactor */

    if (res2 == NULL)
    {
        _arb_hypgeom_legendre_p_ui_zero(t, n, x, xpow, m, K, prec);
        arb_set(res, t);
    }
    else
    {
        _arb_hypgeom_legendre_p_ui_zero(t, n, x, xpow, m, K, prec);
        _arb_hypgeom_legendre_p_ui_zero(u, n - 1, x, xpow, m, K, prec);

        /* P' = n (P[n-1] - x P) / (1 - x^2) */
        arb_submul(u, t, x, prec);
        arb_add_ui(v, v, 1, prec);
        arb_div(u, u, v, prec);
        arb_mul_ui(res2, u, n, prec);
        arb_set(res, t);
    }

    _arb_vec_clear(xpow, m + 1);
    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
}

