/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz_vec.h"
#include "arb_poly.h"
#include "arb_hypgeom.h"

#ifdef __GNUC__
# define pow __builtin_pow
#else
# include <math.h>
#endif

void
arb_hypgeom_rising_ui_jet_rs(arb_ptr res, const arb_t x, ulong n, ulong m, slong len, slong prec)
{
    slong i, j, k, l, m0, xmlen, tlen, ulen, climbs, climbs_max, wp;
    arb_ptr tmp, xpow;
    arb_ptr t, u;
    mp_ptr c;
    TMP_INIT;

    if (len == 0)
        return;

    if (len > n + 1)
    {
        _arb_vec_zero(res + n + 1, len - n - 1);
        len = n + 1;
    }

    if (len == n + 1)
    {
        arb_one(res + n);
        len = n;
    }

    if (n <= 1)
    {
        if (n == 1)
            arb_set_round(res, x, prec);
        return;
    }

    if (len == 1)
    {
        arb_hypgeom_rising_ui_rs(res, x, n, m, prec);
        return;
    }

    TMP_START;

    if (m == 0)
    {
        if (n <= 7)
            m = n;
        else if (n <= 12)
            m = (n + 1) / 2;
        else if (n <= 32)
            m = (n + 2) / 3;
        else
        {
            m0 = n_sqrt(n);
            m = 8 + 0.5 * pow(prec, 0.4);
            m = FLINT_MIN(m, m0);
            m = FLINT_MIN(m, 64);
        }
    }

    wp = ARF_PREC_ADD(prec, FLINT_BIT_COUNT(n));

    climbs_max = FLINT_BIT_COUNT(n - 1) * m;
    c = TMP_ALLOC(sizeof(mp_limb_t) * climbs_max * m);

    /* length of (x+t)^m */
    xmlen = FLINT_MIN(len, m + 1);

    tmp = _arb_vec_init(2 * len + (m + 1) * xmlen);
    t = tmp;
    u = tmp + len;
    xpow = tmp + 2 * len;

    _arb_vec_set_powers(xpow, x, m + 1, wp);

    tlen = 1;

    /* First derivatives */
    for (i = 1; i <= m; i++)
        arb_mul_ui(xpow + (m + 1) + i, xpow + i - 1, i, wp);

    /* Higher derivatives if we need them */
    if (len >= 3)
    {
        fmpz * f = _fmpz_vec_init(len);

        fmpz_one(f + 0);
        fmpz_one(f + 1);

        for (i = 2; i <= m; i++)
        {
            for (j = FLINT_MIN(xmlen - 1, i + 1); j >= 1; j--)
                fmpz_add(f + j, f + j, f + j - 1);

            for (j = 2; j < FLINT_MIN(xmlen, i + 1); j++)
                arb_mul_fmpz(xpow + (m + 1) * j + i, xpow + i - j, f + j, wp);
        }

        _fmpz_vec_clear(f, len);
    }

    for (k = 0; k < n; k += m)
    {
        l = FLINT_MIN(m, n - k);
        climbs = FLINT_BIT_COUNT(k + l - 1) * l;
        climbs = (climbs + FLINT_BITS - 1) / FLINT_BITS;

        ulen = FLINT_MIN(len, l + 1);

        if (l == 1)
        {
            arb_add_ui(u, x, k, wp);
            arb_one(u + 1);
        }
        else
        {
            if (climbs == 1)
            {
                _arb_hypgeom_rising_coeffs_1(c, k, l);
                for (j = 0; j < ulen; j++)
                    arb_dot_ui(u + j, xpow + (m + 1) * j + l, 0, xpow + (m + 1) * j + j, 1, c + j, 1, l - j, wp);
            }
            else if (climbs == 2)
            {
                _arb_hypgeom_rising_coeffs_2(c, k, l);
                for (j = 0; j < ulen; j++)
                    arb_dot_uiui(u + j, xpow + (m + 1) * j + l, 0, xpow + (m + 1) * j + j, 1, c + 2 * j, 1, l - j, wp);
            }
            else
            {
                fmpz * f = (fmpz *) c;

                for (i = 0; i < l; i++)
                    fmpz_init(f + i);

                _arb_hypgeom_rising_coeffs_fmpz(f, k, l);

                for (j = 0; j < ulen; j++)
                    arb_dot_fmpz(u + j, xpow + (m + 1) * j + l, 0, xpow + (m + 1) * j + j, 1, f + j, 1, l - j, wp);

                for (i = 0; i < l; i++)
                    fmpz_clear(f + i);
            }
        }

        if (k == 0)
        {
            tlen = ulen;
            _arb_vec_swap(t, u, ulen);
        }
        else
        {
            _arb_poly_mullow(res, t, tlen, u, ulen, FLINT_MIN(len, tlen + ulen - 1), wp);
            tlen = FLINT_MIN(len, tlen + ulen - 1);
            _arb_vec_swap(t, res, tlen);
        }
    }

    _arb_vec_set_round(res, t, len, prec);

    _arb_vec_clear(tmp, 2 * len + (m + 1) * xmlen);
    TMP_END;
}

