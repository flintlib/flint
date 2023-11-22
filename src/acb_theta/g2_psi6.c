/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_theta.h"

static void
g2_psi6_bits(int * b1, int * b2, int * b3, int * b4, ulong b)
{
    *b4 = b % 2;
    b = b >> 1;
    *b3 = b % 2;
    b = b >> 1;
    *b2 = b % 2;
    b = b >> 1;
    *b1 = b % 2;
}

static slong
g2_psi6_sgn(ulong b, ulong c, ulong d)
{
    slong sgn;
    int b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4;

    g2_psi6_bits(&b1, &b2, &b3, &b4, b);
    g2_psi6_bits(&c1, &c2, &c3, &c4, c);
    g2_psi6_bits(&d1, &d2, &d3, &d4, d);

    sgn = b1 + b2 + c1 + c2 + d1 + d2 + b1 * c1 + b2 * c2 + b4 * c2 + b1 * c3
        - b2 * c4 + b1 * d1 - b3 * d1 + c1 * d1 + b2 * d2 + c2 * d2 + c4 * d2
        + c1 * d3 - b2 * b3 * c1 - b2 * b4 * c2 - b1 * b2 * c3 - b2 * b3 * d1
        - b3 * c1 * d1 - b1 * c3 * d1 - b2 * c3 * d1 - b2 * b4 * d2
        - b4 * c2 * d2 - b1 * b2 * d3 - b1 * c1 * d3 - b2 * c1 * d3;
    sgn = (sgn % 2 == 1 ? -1 : 1);

    return sgn;
}


void
acb_theta_g2_psi6(acb_t res, acb_srcptr th2, slong prec)
{
    slong g = 2;
    ulong ch1, ch2, ch3;
    ulong n = 1 << (2 * g);
    slong sgn;
    acb_t s, t;

    acb_init(s);
    acb_init(t);

    for (ch1 = 0; ch1 < n; ch1++)
    {
        for (ch2 = ch1 + 1; ch2 < n; ch2++)
        {
            for (ch3 = ch2 + 1; ch3 < n; ch3++)
            {
                if (acb_theta_char_is_syzygous(ch1, ch2, ch3, g))
                {
                    sgn = g2_psi6_sgn(ch1, ch2, ch3);
                    acb_mul(t, &th2[ch1], &th2[ch2], prec);
                    acb_mul(t, t, &th2[ch3], prec);
                    acb_sqr(t, t, prec);
                    acb_mul_si(t, t, sgn, prec);
                    acb_add(s, s, t, prec);
                }
            }
        }
    }
    acb_mul_2exp_si(res, s, -2);

    acb_clear(s);
    acb_clear(t);
}
