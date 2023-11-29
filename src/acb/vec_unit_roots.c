/*
    Copyright (C) 2016 Pascal Molin
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
_acb_vec_unit_roots(acb_ptr res, slong n, slong len, slong prec)
{
    int conj = 0;
    slong k, len1, wp;
    acb_t t;

    if (len <= 0)
        return;
    if (n == 0)
    {
        flint_throw(FLINT_ERROR, "\n_acb_vec_unit_roots: need order != 0\n");
    }

    if (n < 0)
    {
        n = -n;
        conj = 1;
    }

    if (n % 4 == 0)
        len1 = FLINT_MIN(len, n / 8 + 1);
    else if (n % 2 == 0)
        len1 = FLINT_MIN(len, n / 4 + 1);
    else
        len1 = FLINT_MIN(len, n / 2 + 1);

    wp = prec + 6 + 2 * FLINT_BIT_COUNT(len1);

    acb_init(t);
    acb_unit_root(t, n, prec);
    _acb_vec_set_powers(res, t, len1, wp);
    acb_clear(t);
    _acb_vec_set_round(res, res, len1, prec);

    if (n % 4 == 0)
    {
        for (k = n / 8 + 1; k <= n / 4 && k < len; k++)
        {
            arb_set(acb_realref(res + k), acb_imagref(res + n / 4 - k));
            arb_set(acb_imagref(res + k), acb_realref(res + n / 4 - k));
        }

        for (k = n / 4 + 1; k <= n / 2 && k < len; k++)
            acb_mul_onei(res + k, res + k - n / 4);
    }
    else if (n % 2 == 0)
    {
        for (k = n / 4 + 1; k <= n / 2 && k < len; k++)
        {
            acb_set(res + k, res + n / 2 - k);
            arb_neg(acb_realref(res + k), acb_realref(res + k));
        }
    }

    for (k = n / 2 + 1; k < len && k < n; k++)
        acb_conj(res + k, res + n - k);

    for (k = n; k < len; k++)
        acb_set(res + k, res - n + k);

    if (conj)
    {
        for (k = 1; k < len; k++)
            acb_conj(res + k, res + k);
    }

}
