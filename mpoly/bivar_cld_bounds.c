/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "long_extras.h"

/*
    for f in F[y][x] of degree n - 1 in x
         on input: l[i] = 1 + deg_y(coeff(f, x^i)) for 0 <= i <= n
        on output: l[i] > deg_y(coeff(f*g_x/g, x^i)) for any divisor g of f

    algo is a 2x scan to find the upper hull and then to fill in the output
*/
void mpoly_bivar_cld_bounds(slong * l, slong n)
{
    slong * P;
    slong Plen = 0;
    slong i, j, x0, y0, x1, y1;
    TMP_INIT;

    FLINT_ASSERT(n > 0);

    TMP_START;

    P = (slong *) TMP_ALLOC(2*n*sizeof(slong));

    n--;

    P[0] = n; P[1] = l[n];
    Plen = 1;
    for (i = n - 1; i >= 0; i--)
    {
        if (l[i] < 1)
            continue;

        x0 = i; y0 = l[i];

        while (Plen >= 2 && !z_mat22_det_is_negative(
                        P[2*(Plen - 1) + 0] - x0, P[2*(Plen - 1) + 1] - y0,
                        P[2*(Plen - 2) + 0] - x0, P[2*(Plen - 2) + 1] - y0))
            Plen--;

        P[2*Plen + 0] = x0;
        P[2*Plen + 1] = y0;
        Plen++;
    }

    i = Plen - 1;
    x0 = P[2*i + 0]; y0 = P[2*i + 1];

    for (j = 1; j <= x0; j++)
        l[j - 1] = j < x0 ? 0 : y0;

    FLINT_ASSERT(j == x0 + 1);

    while (i > 0)
    {
        x1 = P[2*(i - 1) + 0]; y1 = P[2*(i - 1) + 1];
        for ( ; j <= x1; j++)
        {
            ulong t1, t2, t3, t4;
            FLINT_ASSERT(x0 < j && j <= x1);
            umul_ppmm(t1, t2, j - x0, y1);
            umul_ppmm(t3, t4, x1 - j, y0);
            add_ssaaaa(t1, t2, t1, t2, t3, t4);
            udiv_qrnnd(l[j - 1], t3, t1, t2, x1 - x0);
        }
        x0 = x1; y0 = y1;
        i--;
    }

    FLINT_ASSERT(j == n + 1);
    l[j - 1] = 0;

    TMP_END;
}

