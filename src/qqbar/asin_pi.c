/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "qqbar.h"

void best_rational_fast(slong * p, ulong * q, double x, slong N);

int
qqbar_asin_pi(slong * p, ulong * q, const qqbar_t x)
{
    slong deg = qqbar_degree(x);

    *p = 0;
    *q = 0;

    if (deg == 1)
    {
        if (qqbar_is_zero(x))
        {
            *p = 0;
            *q = 1;
            return 1;
        }

        if (qqbar_is_one(x))
        {
            *p = 1;
            *q = 2;
            return 1;
        }

        if (qqbar_is_neg_one(x))
        {
            *p = -1;
            *q = 2;
            return 1;
        }

        if (QQBAR_COEFFS(x)[1] == 2 && QQBAR_COEFFS(x)[0] == -1)
        {
            *p = 1;
            *q = 6;
            return 1;
        }

        if (QQBAR_COEFFS(x)[1] == 2 && QQBAR_COEFFS(x)[0] == 1)
        {
            *p = -1;
            *q = 6;
            return 1;
        }

        return 0;
    }
    else if (deg == 2)
    {
        fmpz a, b, c;

        a = QQBAR_COEFFS(x)[0];
        b = QQBAR_COEFFS(x)[1];
        c = QQBAR_COEFFS(x)[2];

        if (a == -3 && b == 0 && c == 4)
        {
            *p = qqbar_sgn_re(x);
            *q = 3;
            return 1;
        }

        if (a == -1 && b == 0 && c == 2)
        {
            *p = qqbar_sgn_re(x);
            *q = 4;
            return 1;
        }

        if (a == -1 && b == 2 && c == 4)
        {
            *p = (qqbar_sgn_re(x) == 1) ? 1 : -3;
            *q = 10;
            return 1;
        }

        if (a == -1 && b == -2 && c == 4)
        {
            *p = (qqbar_sgn_re(x) == 1) ? 3 : -1;
            *q = 10;
            return 1;
        }

        return 0;
    }
    else if (!qqbar_is_real(x))
    {
        return 0;
    }
    else
    {
        slong degq;
        slong prec;
        arb_t z, pi;
        int res;

        prec = 64;  /* More than enough -- the fractions will only ever be tiny */
        res = 0;

        arb_init(z);
        arb_init(pi);
        qqbar_get_arb(z, x, prec);

        if (arf_cmpabs_2exp_si(arb_midref(z), 0) < 0 && arf_cmpabs_2exp_si(arb_midref(z), -20) > 0)
        {
            arb_asin(z, z, prec);
            arb_const_pi(pi, prec);
            arb_div(z, z, pi, prec);

            best_rational_fast(p, q, arf_get_d(arb_midref(z), ARF_RND_NEAR), 1000000);
            arb_mul_ui(z, z, *q, prec);

            if (arb_contains_si(z, *p))
            {
                if ((*q) % 2 == 1 || (*q) % 4 == 0)
                    degq = n_euler_phi(*q);
                else
                    degq = n_euler_phi(*q) / 2;

                if (deg == degq)
                {
                    qqbar_t v;
                    qqbar_init(v);
                    qqbar_sin_pi(v, *p, *q);
                    res = qqbar_equal(v, x);
                    qqbar_clear(v);
                }
            }
        }

        arb_clear(z);
        arb_clear(pi);

        return res;
    }
}

