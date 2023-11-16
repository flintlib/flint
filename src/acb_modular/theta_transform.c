/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_modular.h"

/* convert theta_{m,n} to theta_i */
static int
swappy1(int m, int n)
{
    m = m & 1;
    n = n & 1;
    if (m == 0 && n == 0) return 2;
    if (m == 0 && n == 1) return 3;
    if (m == 1 && n == 0) return 1;
    return 0;
}

/* extra phase shift picked up: a factor (-1)^n when m is congruent
   to 2,3 (mod 4), and possibly a factor i = sqrt(-1) when converting
   from theta_1 to theta_{1,1} */
static int
swappy2(int m, int n)
{
    m = m & 3;
    n = n & 1;
    if (m == 1 && n == 1) return 2; /* i */
    if (m == 2 && n == 1) return 4; /* (-1)^n */
    if (m == 3 && n == 1) return 6; /* (-1)^n * i */
    return 0;
}

void
acb_modular_theta_transform(int * R, int * S, int * C, const psl2z_t g)
{
    R[0] = 0;
    R[1] = 0;
    R[2] = 0;
    R[3] = 0;

    S[0] = 0;
    S[1] = 1;
    S[2] = 2;
    S[3] = 3;

    if (fmpz_is_zero(&g->c))
    {
        C[0] = 0;

        if (fmpz_is_odd(&g->b))
        {
            S[2] = 3;
            S[3] = 2;
        }

        /* -b mod 8 */
        R[0] = (- (int) fmpz_fdiv_ui(&g->b, 8)) & 7;
        R[1] = R[0];
    }
    else
    {
        int a, b, c, d, e1, e2;
        psl2z_t h;

        psl2z_init(h);
        psl2z_inv(h, g);

        e1 = acb_modular_epsilon_arg(h);
        e2 = acb_modular_epsilon_arg(g);

        psl2z_clear(h);

        C[0] = 1;

        a = fmpz_fdiv_ui(&g->a, 8);
        b = fmpz_fdiv_ui(&g->b, 8);
        c = fmpz_fdiv_ui(&g->c, 8);
        d = fmpz_fdiv_ui(&g->d, 8);

        R[0] = e1 + 1;
        R[1] = -e2 + 5 + (2 - c) * a;
        R[2] = -e2 + 4 + (c - d - 2) * (b - a);
        R[3] = -e2 + 3 - (2 + d) * b;

        S[1] = swappy1(1 - c, 1 + a);
        R[1] += swappy2(1 - c, 1 + a);

        S[2] = swappy1(1 + d - c, 1 - b + a);
        R[2] += swappy2(1 + d - c, 1 - b + a);

        S[3] = swappy1(1 + d, 1 - b);
        R[3] += swappy2(1 + d, 1 - b);

        /* floor mod by 8 */
        R[0] &= 7;
        R[1] &= 7;
        R[2] &= 7;
        R[3] &= 7;
    }
}

