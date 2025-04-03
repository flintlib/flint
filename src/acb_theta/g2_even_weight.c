/*
    Copyright (C) 2025 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_theta.h"

static int
acb_theta_char_is_goepel(ulong ch1, ulong ch2, ulong ch3, ulong ch4, slong g)
{
    if (ch1 == ch2 || ch1 == ch3 || ch1 == ch4
        || ch2 == ch3 || ch2 == ch4 || ch3 == ch4)
    {
        return 0;
    }

    return acb_theta_char_is_even(ch1, g)
        && acb_theta_char_is_even(ch2, g)
        && acb_theta_char_is_even(ch3, g)
        && acb_theta_char_is_even(ch4, g)
        && ((ch1 ^ ch2 ^ ch3 ^ ch4) == 0);
}

static int
acb_theta_char_is_syzygous(ulong ch1, ulong ch2, ulong ch3, slong g)
{
    return acb_theta_char_is_goepel(ch1, ch2, ch3, ch1 ^ ch2 ^ ch3, g);
}

static void
g2_psi6_bits(int * b1, int * b2, int * b3, int * b4, ulong b)
{
    *b1 = acb_theta_char_bit(b, 0, 4);
    *b2 = acb_theta_char_bit(b, 1, 4);
    *b3 = acb_theta_char_bit(b, 2, 4);
    *b4 = acb_theta_char_bit(b, 3, 4);
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

static void
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

static void
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

void
acb_theta_g2_even_weight(acb_t psi4, acb_t psi6, acb_t chi10, acb_t chi12,
    acb_srcptr th2, slong prec)
{
    slong g = 2;
    slong n = 1 << (2 * g);
    ulong ab;
    acb_t s, t;

    acb_init(s);
    acb_init(t);

    for (ab = 0; ab < (1 << (2 * g)); ab++)
    {
        if (acb_theta_char_is_even(ab, g))
        {
            acb_pow_ui(t, &th2[ab], 4, prec);
            acb_add(s, s, t, prec);
        }
    }
    acb_mul_2exp_si(psi4, s, -2);

    acb_one(t);
    for (ab = 0; ab < n; ab++)
    {
        if (acb_theta_char_is_even(ab, g))
        {
            acb_mul(t, t, &th2[ab], prec);
        }
    }
    acb_mul_2exp_si(chi10, t, -12);

    acb_theta_g2_psi6(psi6, th2, prec);
    acb_theta_g2_chi12(chi12, th2, prec);

    acb_clear(t);
    acb_clear(s);
}
