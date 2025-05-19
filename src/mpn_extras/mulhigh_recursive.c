/*
    Copyright (C) 2024 Albin Ahlbäck
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"

mp_limb_t
_flint_mpn_mulhigh_n_recursive(mp_ptr r, mp_srcptr x, mp_srcptr y, mp_size_t n)
{
    if (FLINT_HAVE_MULHIGH_FUNC(n))
    {
        return flint_mpn_mulhigh_func_tab[n](r, x, y);
    }
    else if (n <= 2 * FLINT_MPN_MULHIGH_BEST_TAB_N)
    {
        mp_limb_t t[2 * FLINT_MPN_MULHIGH_BEST_TAB_N];

        mp_size_t m1 = n - FLINT_MPN_MULHIGH_BEST_TAB_N;
        mp_size_t m2 = FLINT_MPN_MULHIGH_BEST_TAB_N;
        mp_limb_t cy, lo, w0, w1, w2;

        FLINT_ASSERT(FLINT_MPN_MULHIGH_BEST_TAB_N <= FLINT_MPN_MULHIGH_FUNC_TAB_WIDTH);

        flint_mpn_mul(r, x + m1, m2, y + m2, m1);
        w0 = flint_mpn_mulhigh_n(t, x + m1, y, m2);
        cy = mpn_add_n(r, r, t, m2);
        MPN_INCR_U(r + m2, m1, cy);
        w1 = flint_mpn_mulhigh_n(t, x, y + m2, m1);
        cy = mpn_add_n(r, r, t, m1);
        MPN_INCR_U(r + m1, m2, cy);

        umul_ppmm(w2, lo, x[m1 - 1], y[m2 - 1]);
        add_ssaaaa(cy, w0, 0, w0, 0, w1);
        add_ssaaaa(cy, w0, cy, w0, 0, w2);
        MPN_INCR_U(r, n, cy);

        return w0;
    }
    else
    {
        mp_ptr t;
        mp_limb_t cy, lo, w0, w1, w2;
        mp_size_t m1 = n / 2;
        mp_size_t m2 = n - m1;
        TMP_INIT;

        TMP_START;
        t = TMP_ALLOC(sizeof(mp_limb_t) * n);

        if (n % 2 == 0)
        {
            flint_mpn_mul_n(r, x + m1, y + m1, m1);
            w0 = _flint_mpn_mulhigh_n_recursive(t, x + m1, y, m1);
            cy = mpn_add_n(r, r, t, m1);
            w1 = _flint_mpn_mulhigh_n_recursive(t, x, y + m1, m1);
            cy += mpn_add_n(r, r, t, m1);
            MPN_INCR_U(r + m1, m1, cy);
        }
        else
        {
            flint_mpn_mul(r, x + m1, m2, y + m2, m1);
            w0 = _flint_mpn_mulhigh_n_recursive(t, x + m1, y, m2);
            cy = mpn_add_n(r, r, t, m2);
            MPN_INCR_U(r + m2, m1, cy);
            w1 = _flint_mpn_mulhigh_n_recursive(t, x, y + m2, m1);
            cy = mpn_add_n(r, r, t, m1);
            MPN_INCR_U(r + m1, m2, cy);
        }

        umul_ppmm(w2, lo, x[m1 - 1], y[m2 - 1]);
        add_ssaaaa(cy, w0, 0, w0, 0, w1);
        add_ssaaaa(cy, w0, cy, w0, 0, w2);
        MPN_INCR_U(r, n, cy);

        TMP_END;

        return w0;
    }
}
