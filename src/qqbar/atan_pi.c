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

#ifdef __GNUC__
# define fabs __builtin_fabs
# define floor __builtin_floor
#else
# include <math.h>
#endif

void
best_rational_fast(slong * p, ulong * q, double x, slong N)
{
    slong a, b, c, d;
    double m, t, u, eps;

    if (x > 1.0 || x < 0.0)
    {
        double n = floor(x);

        best_rational_fast(p, q, x - n, N);
        *p = *p + n * (*q);
        return;
    }

    a = 0; b = 1;
    c = 1; d = 1;

    eps = 0.1 / N;

    if (fabs(x) < eps)
    {
        *p = 0;
        *q = 1;
        return;
    }

    while (b <= N && d <= N)
    {
        m = (a + c) / ((double) (b + d));

        if (fabs(m - x) < eps)
        {
            if (b + d <= N)
            {
                *p = a + c;
                *q = b + d;
            }
            else if (d > b)
            {
                *p = c;
                *q = d;
            }
            else
            {
                *p = a;
                *q = b;
            }

            return;
        }
        else
        {
            t = a + c;
            u = b + d;

            if (x > m)
            {
                a = t;
                b = u;
            }
            else
            {
                c = t;
                d = u;
            }
        }
    }

    if (b > N)
    {
        *p = c;
        *q = d;
    }
    else
    {
        *p = a;
        *q = b;
    }
}

int
qqbar_atan_pi(slong * p, ulong * q, const qqbar_t x)
{
    slong deg = qqbar_degree(x);

    *p = 0;
    *q = 1;

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
            *q = 4;
            return 1;
        }

        if (qqbar_is_neg_one(x))
        {
            *p = -1;
            *q = 4;
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

        if (a == -3 && b == 0 && c == 1)
        {
            *p = qqbar_sgn_re(x);
            *q = 3;
            return 1;
        }

        if (a == -1 && b == 0 && c == 3)
        {
            *p = qqbar_sgn_re(x);
            *q = 6;
            return 1;
        }

        if (a == -1 && b == 2 && c == 1)
        {
            *p = (qqbar_sgn_re(x) == 1) ? 1 : -3;
            *q = 8;
            return 1;
        }

        if (a == -1 && b == -2 && c == 1)
        {
            *p = (qqbar_sgn_re(x) == 1) ? 3 : -1;
            *q = 8;
            return 1;
        }

        if (a == 1 && b == -4 && c == 1)
        {
            /* root is ~0.267 or ~3.73 -- accuracy should not be that bad */
            if (arb_contains_si(acb_realref(QQBAR_ENCLOSURE(x)), 1))
                flint_throw(FLINT_ERROR, "(%s)\n", __func__);

           *p = (arf_cmpabs_2exp_si(arb_midref(acb_realref(QQBAR_ENCLOSURE(x))), 0) < 0) ? 1 : 5;
            *q = 12;
            return 1;
        }

        if (a == 1 && b == 4 && c == 1)
        {
            if (arb_contains_si(acb_realref(QQBAR_ENCLOSURE(x)), -1))
                flint_throw(FLINT_ERROR, "(%s)\n", __func__);
            *p = (arf_cmpabs_2exp_si(arb_midref(acb_realref(QQBAR_ENCLOSURE(x))), 0) < 0) ? -1 : -5;
            *q = 12;
            return 1;
        }

        return 0;
    }
    else if ((deg % 2 != 0) || !qqbar_is_real(x))
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

        if (arf_cmpabs_2exp_si(arb_midref(z), 20) < 0 && arf_cmpabs_2exp_si(arb_midref(z), -20) > 0)
        {
            arb_atan(z, z, prec);
            arb_const_pi(pi, prec);
            arb_div(z, z, pi, prec);

            best_rational_fast(p, q, arf_get_d(arb_midref(z), ARF_RND_NEAR), 1000000);
            arb_mul_ui(z, z, *q, prec);

            if (arb_contains_si(z, *p))
            {
                if ((*q) % 4 == 0)
                    degq = n_euler_phi(*q) / 2;
                else
                    degq = n_euler_phi(*q);

                if (deg == degq)
                {
                    qqbar_t v;
                    qqbar_init(v);
                    qqbar_tan_pi(v, *p, *q);
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

