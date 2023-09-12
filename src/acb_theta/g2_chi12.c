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
acb_theta_g2_chi12(acb_t r, acb_srcptr th2, slong prec)
{
    slong g = 2;
    ulong ch1, ch2, ch3, ch4, ab;
    ulong n = 1 << (2 * g);
    acb_t res, aux;

    acb_init(res);
    acb_init(aux);

    for (ch1 = 0; ch1 < n; ch1++)
    {
        for (ch2 = ch1 + 1; ch2 < n; ch2++)
        {
            for (ch3 = ch2 + 1; ch3 < n; ch3++)
            {
                for (ch4 = ch3 + 1; ch4 < n; ch4++)
                {
                    if (acb_theta_char_is_goepel(ch1, ch2, ch3, ch4, g))
                    {
                        acb_one(aux);
                        for (ab = 0; ab < n; ab++)
                        {
                            if (acb_theta_char_is_even(ab, g)
                                && (ab != ch1)
                                && (ab != ch2)
                                && (ab != ch3)
                                && (ab != ch4))
                            {
                                acb_mul(aux, aux, &th2[ab], prec);
                            }
                        }
                        acb_sqr(aux, aux, prec);
                        acb_add(res, res, aux, prec);
                    }
                }
            }
        }
    }
    acb_mul_2exp_si(r, res, -15);

    acb_clear(res);
    acb_clear(aux);
}
