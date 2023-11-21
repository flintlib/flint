/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"

#ifdef __GNUC__
# define pow __builtin_pow
#else
# include <math.h>
#endif

void
acb_hypgeom_rising_ui_rs(acb_t res, const acb_t x, ulong n, ulong m, slong prec)
{
    slong i, k, l, m0, climbs, climbs_max, wp;
    acb_ptr xpow;
    acb_t t, u;
    mp_ptr c;
    TMP_INIT;

    if (n <= 1)
    {
        if (n == 0)
            acb_one(res);
        else
            acb_set_round(res, x, prec);
        return;
    }

    TMP_START;

    if (m == 0 || m == -1)
    {
        if (n <= 6)
            m = 2;
        else if (n <= 16)
            m = 4;
        else if (n <= 40)
            m = 6;
        else
        {
            m0 = n_sqrt(n);
            m = 8 + 0.27 * pow(1.5 * FLINT_MAX(0, prec - 1024), 0.4);
            m = FLINT_MIN(m, m0);
            m = FLINT_MIN(m, 64);
        }
    }

    wp = ARF_PREC_ADD(prec, FLINT_BIT_COUNT(n));

    climbs_max = FLINT_BIT_COUNT(n - 1) * m;
    c = TMP_ALLOC(sizeof(mp_limb_t) * climbs_max * m);

    xpow = _acb_vec_init(m + 1);
    _acb_vec_set_powers(xpow, x, m + 1, wp);
    acb_init(t);
    acb_init(u);

    for (k = 0; k < n; k += m)
    {
        l = FLINT_MIN(m, n - k);
        climbs = FLINT_BIT_COUNT(k + l - 1) * l;
        climbs = (climbs + FLINT_BITS - 1) / FLINT_BITS;

        /* assumes l >= 2 */
        if (l == 1)
        {
            acb_add_ui(u, x, k, wp);
        }
        else
        {
            if (climbs == 1)
            {
                _arb_hypgeom_rising_coeffs_1(c, k, l);
                acb_dot_ui(u, xpow + l, 0, xpow, 1, c, 1, l, wp);
            }
            else if (climbs == 2)
            {
                _arb_hypgeom_rising_coeffs_2(c, k, l);
                acb_dot_uiui(u, xpow + l, 0, xpow, 1, c, 1, l, wp);
            }
            else
            {
                fmpz * f = (fmpz *) c;

                for (i = 0; i < l; i++)
                    fmpz_init(f + i);

                _arb_hypgeom_rising_coeffs_fmpz(f, k, l);

                acb_dot_fmpz(u, xpow + l, 0, xpow, 1, f, 1, l, wp);

                for (i = 0; i < l; i++)
                    fmpz_clear(f + i);
            }
        }

        if (k == 0)
            acb_swap(t, u);
        else
            acb_mul(t, t, u, wp);
    }

    acb_set_round(res, t, prec);

    acb_clear(t);
    acb_clear(u);
    _acb_vec_clear(xpow, m + 1);
    TMP_END;
}

