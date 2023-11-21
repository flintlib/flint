/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mag.h"

void
mag_bin_uiui(mag_t res, ulong n, ulong k)
{
    mag_t t;

    if (k > n)
    {
        mag_zero(res);
        return;
    }

    if (k > n / 2)
        k = n - k;

    if (k == 0)
    {
        mag_one(res);
        return;
    }

    if (k == 1)
    {
        mag_set_ui(res, n);
        return;
    }

    mag_init(t);

    if (n < 256 && k < 256)
    {
        /* using accurate (lookup table) factorials */
        mag_fac_ui(res, n);
        mag_rfac_ui(t, k);
        mag_mul(res, res, t);
        mag_rfac_ui(t, n - k);
        mag_mul(res, res, t);
    }
    else
    {
        /* binary entropy bound (n/(n-k))^(n-k) (n/k)^k = n^n / (k^k (n-k)^(n-k)) */
        mag_set_ui(res, n);
        mag_div_ui(res, res, n - k);
        mag_pow_ui(res, res, n - k);
        /* (n/(n-k))^(n-k) is also bounded by exp(k) */
        mag_set_ui(t, k);
        mag_exp(t, t);
        mag_min(res, res, t);
        mag_set_ui(t, n);
        mag_div_ui(t, t, k);
        mag_pow_ui(t, t, k);
        mag_mul(res, res, t);
    }

    mag_clear(t);
}

