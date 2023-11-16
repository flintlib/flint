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
sum_rs_inner(arb_t s, arb_srcptr xpow, slong m, ulong n, slong K, ulong prime, slong prec)
{
    slong j, k, khi, klo, u, r;
    fmpz * c;

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
            ulong aa = (n - j + 1 - prime);
            ulong bb = (n + j + prime);

            if (j == klo)
                fmpz_ui_mul_ui(c + khi - j, aa, bb);
            else
                fmpz_mul2_uiui(c + khi - j, c + khi - j + 1, aa, bb);
        }

        for (j = khi; j >= klo; j--)
        {
            ulong aa = (j);
            ulong bb = (j + prime);

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

void
arb_hypgeom_legendre_p_ui_one(arb_t res, arb_t res_prime, ulong n,
                        const arb_t x, slong K, slong prec)
{
    arb_t s, v;
    arb_ptr xpow;
    slong m;
    mag_t u, a, t, err;

    if (n == 0)
    {
        if (res != NULL) arb_one(res);
        if (res_prime != NULL) arb_zero(res_prime);
        return;
    }

    /* overflow protection */
    if (n > UWORD_MAX / 4)
    {
        if (res != NULL) arb_indeterminate(res);
        if (res_prime != NULL) arb_indeterminate(res_prime);
    }

    arb_init(v);
    arb_init(s);
    mag_init(u);
    mag_init(a);
    mag_init(t);
    mag_init(err);

    K = FLINT_MIN(K, n + 1);

    if (res != NULL && res_prime != NULL)
        m = n_sqrt(2 * K);
    else
        m = n_sqrt(K);

    xpow = _arb_vec_init(m + 1);

    arb_sub_ui(v, x, 1, prec);
    arb_mul_2exp_si(v, v, -1);
    _arb_vec_set_powers(xpow, v, m + 1, prec);

    /* truncating */
    if (K < n + 1)
    {
        arb_get_mag(u, v);
        mag_mul_ui(t, u, n - K);
        mag_mul_ui(t, t, n + K + 1);
        mag_div_ui(t, t, K + 1);
        mag_div_ui(t, t, K + 1);
        mag_geom_series(t, t, 0);
        mag_pow_ui(u, u, K);
        mag_mul(u, u, t);
    }

    if (res != NULL)
    {
        sum_rs_inner(s, xpow, m, n, K, 0, prec);
        arb_add_ui(res, s, 1, prec);

        if (K < n + 1)
        {
            mag_set(err, u);
            mag_bin_uiui(t, n, K);
            mag_mul(err, err, t);
            mag_bin_uiui(t, n + K, K);
            mag_mul(err, err, t);
            arb_add_error_mag(res, err);
        }
    }

    if (res_prime != NULL)
    {
        K = FLINT_MIN(K, n);
        sum_rs_inner(s, xpow, m, n, K, 1, prec);
        arb_add_ui(res_prime, s, 1, prec);
        arb_mul_ui(res_prime, res_prime, n, prec);
        arb_mul_ui(res_prime, res_prime, n + 1, prec);
        arb_mul_2exp_si(res_prime, res_prime, -1);

        /* truncating */
        if (K < n)
        {
            mag_set(err, u);
            mag_bin_uiui(t, n, K + 1);
            mag_mul(err, err, t);
            mag_bin_uiui(t, n + K + 1, K + 1);
            mag_mul(err, err, t);
            mag_mul_ui(err, err, n);
            arb_add_error_mag(res_prime, err);
        }

    }

    _arb_vec_clear(xpow, m + 1);
    arb_clear(s);
    arb_clear(v);
    mag_clear(u);
    mag_clear(a);
    mag_clear(t);
    mag_clear(err);
}

