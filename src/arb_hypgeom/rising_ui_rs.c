/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "arb_hypgeom.h"

#ifdef __GNUC__
# define pow __builtin_pow
#else
# include <math.h>
#endif

void
_arb_hypgeom_rising_coeffs_1(ulong * c, ulong k, slong l)
{
    slong i, j;
    ulong d;

    if (l < 2)
        flint_throw(FLINT_ERROR, "(%s): l < 2\n", __func__);

    c[0] = k * (k + 1);
    c[1] = 2 * k + 1;

    for (i = 2; i < l; i++)
    {
        d = k + i;

        c[i] = c[i - 1] + d;

        for (j = i - 1; j >= 1; j--)
            c[j] = c[j] * d + c[j - 1];

        c[0] = c[0] * d;
    }
}

void
_arb_hypgeom_rising_coeffs_2(ulong * c, ulong k, slong l)
{
    slong i, j;
    ulong d;
    mp_limb_t hi, lo;

    if (l < 2)
        flint_throw(FLINT_ERROR, "(%s): l < 2\n", __func__);

    umul_ppmm(c[1], c[0], k, k + 1);

    c[2] = 2 * k + 1;
    c[3] = 0;

    for (i = 2; i < l; i++)
    {
        d = k + i;

        add_ssaaaa(c[2 * i + 1], c[2 * i], c[2 * i - 1], c[2 * i - 2], 0, d);

        for (j = i - 1; j >= 1; j--)
        {
            umul_ppmm(hi, lo, c[2 * j], d);
            hi += c[2 * j + 1] * d;
            add_ssaaaa(c[2 * j + 1], c[2 * j], hi, lo, c[2 * j - 1], c[2 * j - 2]);
        }

        umul_ppmm(hi, lo, c[0], d);
        c[0] = lo;
        c[1] = c[1] * d + hi;
    }
}

void
_arb_hypgeom_rising_coeffs_fmpz(fmpz * c, ulong k, slong l)
{
    slong i, j;
    ulong d;

    if (l < 2)
        flint_throw(FLINT_ERROR, "(%s): l < 2\n", __func__);

    fmpz_set_ui(c + 0, k);
    fmpz_mul_ui(c + 0, c + 0, k + 1);
    fmpz_set_ui(c + 1, 2 * k + 1);

    for (i = 2; i < l; i++)
    {
        d = k + i;

        fmpz_add_ui(c + i, c + i - 1, d);

        for (j = i - 1; j >= 1; j--)
        {
            fmpz_mul_ui(c + j, c + j, d);
            fmpz_add(c + j, c + j, c + j - 1);
        }

        fmpz_mul_ui(c + 0, c + 0, d);
    }
}

void
arb_hypgeom_rising_ui_rs(arb_t res, const arb_t x, ulong n, ulong m, slong prec)
{
    slong i, k, l, m0, climbs, climbs_max, wp;
    arb_ptr xpow;
    arb_t t, u;
    mp_ptr c;
    TMP_INIT;

    if (n <= 1)
    {
        if (n == 0)
            arb_one(res);
        else
            arb_set_round(res, x, prec);
        return;
    }

    TMP_START;

    if (m == 0)
    {
        if (n <= 6)
            m = 1 + (prec >= 1024);
        else if (n <= 16)
            m = 4;
        else if (n <= 50)
            m = 6;
        else
        {
            m0 = n_sqrt(n);
            m = 8 + 0.2 * pow(FLINT_MAX(0, prec - 4096), 0.4);
            m = FLINT_MIN(m, m0);
            m = FLINT_MIN(m, 60);
        }
    }

    wp = ARF_PREC_ADD(prec, FLINT_BIT_COUNT(n));

    climbs_max = FLINT_BIT_COUNT(n - 1) * m;
    c = TMP_ALLOC(sizeof(mp_limb_t) * climbs_max * m);

    xpow = _arb_vec_init(m + 1);
    _arb_vec_set_powers(xpow, x, m + 1, wp);
    arb_init(t);
    arb_init(u);

    for (k = 0; k < n; k += m)
    {
        l = FLINT_MIN(m, n - k);
        climbs = FLINT_BIT_COUNT(k + l - 1) * l;
        climbs = (climbs + FLINT_BITS - 1) / FLINT_BITS;

        /* assumes l >= 2 */
        if (l == 1)
        {
            arb_add_ui(u, x, k, wp);
        }
        else
        {
            if (climbs == 1)
            {
                _arb_hypgeom_rising_coeffs_1(c, k, l);
                arb_dot_ui(u, xpow + l, 0, xpow, 1, c, 1, l, wp);
            }
            else if (climbs == 2)
            {
                _arb_hypgeom_rising_coeffs_2(c, k, l);
                arb_dot_uiui(u, xpow + l, 0, xpow, 1, c, 1, l, wp);
            }
            else
            {
                fmpz * f = (fmpz *) c;

                for (i = 0; i < l; i++)
                    fmpz_init(f + i);

                _arb_hypgeom_rising_coeffs_fmpz(f, k, l);

                arb_dot_fmpz(u, xpow + l, 0, xpow, 1, f, 1, l, wp);

                for (i = 0; i < l; i++)
                    fmpz_clear(f + i);
            }
        }

        if (k == 0)
            arb_swap(t, u);
        else
            arb_mul(t, t, u, wp);
    }

    arb_set_round(res, t, prec);

    arb_clear(t);
    arb_clear(u);
    _arb_vec_clear(xpow, m + 1);
    TMP_END;
}

