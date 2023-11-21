/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_modular.h"

#ifdef __GNUC__
# define floor __builtin_floor
#else
# include <math.h>
#endif

static inline int
d_is_ok(double x)
{
    return (x > -1e15) && (x < 1e15);
}

void
acb_modular_fundamental_domain_approx_d(psl2z_t g,
    double x, double y, double one_minus_eps)
{
    double a, b, c, d, aa, bb, t;
    int i;

    a = d = 1;
    b = c = 0;

    for (i = 0; i < 20; i++)
    {
        if (!d_is_ok(x) || !d_is_ok(y) || !(y > 0.0))
        {
            psl2z_one(g);
            return;
        }

        /* shift */
        if (x < -0.5 || x > 0.5)
        {
            t = floor(x + 0.5);
            x -= t;
            a -= t * c;
            b -= t * d;

            /* too large to guarantee exactness */
            if (!d_is_ok(a) || !d_is_ok(b))
            {
                psl2z_one(g);
                return;
            }

            continue;
        }

        t = x*x + y*y;

        /* can't divide by a too small number */
        if (t < 1e-30)
        {
            psl2z_one(g);
            break;
        }

        /* inversion */
        if (t < one_minus_eps)
        {
            t = 1.0 / t;
            x *= -t;
            y *= t;
            aa = a;
            bb = b;
            a = -c;
            b = -d;
            c = aa;
            d = bb;
            continue;
        }

        /* we're good */
        break;
    }

    if (c < 0 || (c == 0 && d < 0))
    {
        a = -a;
        b = -b;
        c = -c;
        d = -d;
    }

    fmpz_set_d(&g->a, a);
    fmpz_set_d(&g->b, b);
    fmpz_set_d(&g->c, c);
    fmpz_set_d(&g->d, d);
}

