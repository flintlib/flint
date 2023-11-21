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

void
acb_theta_g2_chi12(acb_t res, acb_srcptr th2, slong prec)
{
    slong g = 2;
    ulong ch1, ch2, ch3, ch4, ab;
    ulong n = 1 << (2 * g);
    acb_t s, t;

    acb_init(s);
    acb_init(t);

    for (ch1 = 0; ch1 < n; ch1++)
    {
        for (ch2 = ch1 + 1; ch2 < n; ch2++)
        {
            for (ch3 = ch2 + 1; ch3 < n; ch3++)
            {
                ch4 = ch1 ^ ch2 ^ ch3;
                if (acb_theta_char_is_goepel(ch1, ch2, ch3, ch4, g))
                {
                    acb_one(t);
                    for (ab = 0; ab < n; ab++)
                    {
                        if (acb_theta_char_is_even(ab, g)
                            && (ab != ch1)
                            && (ab != ch2)
                            && (ab != ch3)
                            && (ab != ch4))
                        {
                            acb_mul(t, t, &th2[ab], prec);
                        }
                    }
                    acb_sqr(t, t, prec);
                    acb_add(s, s, t, prec);
                }
            }
        }
    }
    acb_mul_2exp_si(res, s, -15);

    acb_clear(s);
    acb_clear(t);
}
