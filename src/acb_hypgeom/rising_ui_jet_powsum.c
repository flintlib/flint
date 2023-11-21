/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "fmpz_vec.h"
#include "acb.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_rising_ui_jet_powsum(acb_ptr res, const acb_t x, ulong n, slong len, slong prec)
{
    slong i, j, k, wp;
    acb_ptr xpow;
    TMP_INIT;

    if (len == 0)
        return;

    if (len > n + 1)
    {
        _acb_vec_zero(res + n + 1, len - n - 1);
        len = n + 1;
    }

    if (len == n + 1)
    {
        acb_one(res + n);
        len = n;
    }

    if (n <= 1)
    {
        if (n == 1)
            acb_set_round(res, x, prec);
        return;
    }

    if (len == 1)
    {
        acb_hypgeom_rising_ui_rs(res, x, n, 0, prec);
        return;
    }

    if (n == 2)
    {
        acb_mul_2exp_si(res + 1, x, 1);
        acb_add_ui(res + 1, res + 1, 1, prec);
        acb_mul(res, x, x, prec + 4);
        acb_add(res, res, x, prec);
        return;
    }

    if (n <= 12 || (FLINT_BITS == 64 && n <= 20))
    {
        mp_ptr c;
        TMP_START;

        wp = ARF_PREC_ADD(prec, FLINT_BIT_COUNT(n));
        c = TMP_ALLOC(sizeof(mp_limb_t) * (n + 1) * len);

        _nmod_vec_zero(c, (n + 1) * len);

        c[0] = 0;
        c[1] = 1;
        c[(n + 1) + 0] = 1;

        for (i = 2; i <= n; i++)
        {
            for (j = FLINT_MIN(len - 1, i); j >= 0; j--)
            {
                slong ln, pos;

                ln = i + 1 - j;
                pos = (n + 1) * j;
                if (i == j)
                {
                    c[pos] = 1;
                }
                else
                {
                    c[pos + ln - 1] = c[pos + ln - 2];
                    for (k = ln - 2; k >= 1; k--)
                        c[pos + k] = c[pos + k] * (i - 1) + c[pos + k - 1];
                    c[pos + 0] *= (i - 1);
                    if (j != 0)
                        for (k = ln - 1; k >= 0; k--)
                            c[pos + k] += c[pos - (n + 1) + k];
                }
            }
        }

        xpow = _acb_vec_init(n + 1);
        _acb_vec_set_powers(xpow, x, n + 1, wp);

        acb_dot_ui(res, NULL, 0, xpow + 1, 1, c + 1, 1, n, prec);

        for (i = 1; i < len; i++)
            acb_dot_ui(res + i, NULL, 0, xpow, 1, c + (n + 1) * i, 1, n + 1 - i, prec);

        _acb_vec_clear(xpow, n + 1);
        TMP_END;
    }
    else
    {
        fmpz * c;

        wp = ARF_PREC_ADD(prec, FLINT_BIT_COUNT(n));
        c = _fmpz_vec_init((n + 1) * len);

        fmpz_one(c + 1);
        fmpz_one(c + (n + 1) + 0);

        for (i = 2; i <= n; i++)
        {
            for (j = FLINT_MIN(len - 1, i); j >= 0; j--)
            {
                slong ln, pos;

                ln = i + 1 - j;
                pos = (n + 1) * j;
                if (i == j)
                {
                    fmpz_one(c + pos);
                }
                else
                {
                    fmpz_set(c + pos + ln - 1, c + pos + ln - 2);
                    for (k = ln - 2; k >= 1; k--)
                    {
                        fmpz_mul_ui(c + pos + k, c + pos + k, i - 1);
                        fmpz_add(c + pos + k, c + pos + k, c + pos + k - 1);
                    }

                    fmpz_mul_ui(c + pos + 0, c + pos + 0, i - 1);
                    if (j != 0)
                        for (k = ln - 1; k >= 0; k--)
                            fmpz_add(c + pos + k, c + pos + k, c + pos - (n + 1) + k);
                }
            }
        }

        xpow = _acb_vec_init(n + 1);
        _acb_vec_set_powers(xpow, x, n + 1, wp);

        acb_dot_fmpz(res, NULL, 0, xpow + 1, 1, c + 1, 1, n, prec);

        for (i = 1; i < len; i++)
            acb_dot_fmpz(res + i, NULL, 0, xpow, 1, c + (n + 1) * i, 1, n + 1 - i, prec);

        _acb_vec_clear(xpow, n + 1);
        _fmpz_vec_clear(c, (n + 1) * len);
    }
}

