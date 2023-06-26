/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_dupl_z(acb_ptr r, acb_srcptr th, slong g, slong prec)
{
    slong n = 1 << g;
    acb_ptr v;
    acb_t c;
    ulong a, b, aa, bb;
    slong s;

    v = _acb_vec_init(n * n);
    acb_init(c);

    for (a = 0; a < n; a++)
    {
        for (b = 0; b < n; b++)
        {
            for (aa = 0; aa < n; aa++)
            {
                for (bb = 0; bb < n; bb++)
                {
                    acb_mul(c, &th[n * aa + bb], &th[(n * aa + bb) ^ b], prec);
                    acb_mul(c, c, &th[(n * aa + bb) ^ (n * a)], prec);
                    acb_mul(c, c, &th[(n * aa + bb) ^ (n * a + b)], prec);

                    s = (acb_theta_char_dot(a, bb, g)
                         + acb_theta_char_dot(a, bb & b, g)) % 2;
                    if (s == 1)
                    {
                        acb_neg(c, c);
                    }
                    acb_add(&v[n * a + b], &v[n * a + b], c, prec);
                }
            }
            acb_div(&v[n * a + b], &v[n * a + b], &th[n * n + n * a], prec);
            acb_div(&v[n * a + b], &v[n * a + b], &th[n * n + b], prec);
        }
    }
    _acb_vec_scalar_mul_2exp_si(v, v, n * n, -g);
    _acb_vec_scalar_div(v, v, n * n, &th[n * n], prec);

    _acb_vec_set(r, v, n * n);
    _acb_vec_set(r + n * n, th + n * n, n * n);

    _acb_vec_clear(v, n * n);
    acb_clear(c);
}
