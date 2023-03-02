/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "d_vec.h"

double
_d_vec_dot_thrice(const double *vec1, const double *vec2, slong len2,
                  double *err)
{
    int i, j;
    double p, h, a1, a2, b1, b2, c, res = 0, g;
    double *r;
    ulong factor = (UWORD(1) << 27) + 1;

    if (len2 == 0)
    {
        *err = 0;
        return 0;
    }

    r = _d_vec_init(2 * len2);

    p = vec1[0] * vec2[0];
    c = factor * vec1[0];
    a1 = c - (c - vec1[0]);
    a2 = vec1[0] - a1;
    c = factor * vec2[0];
    b1 = c - (c - vec2[0]);
    b2 = vec2[0] - b1;
    r[0] = a2 * b2 - (((p - a1 * b1) - a2 * b1) - a1 * b2);

    for (i = 1; i < len2; i++)
    {
        h = vec1[i] * vec2[i];
        c = factor * vec1[i];
        a1 = c - (c - vec1[i]);
        a2 = vec1[i] - a1;
        c = factor * vec2[i];
        b1 = c - (c - vec2[i]);
        b2 = vec2[i] - b1;
        r[i] = a2 * b2 - (((h - a1 * b1) - a2 * b1) - a1 * b2);

        a1 = p;
        p = p + h;
        c = p - a1;
        r[len2 + i - 1] = (a1 - (p - c)) + (h - c);
    }

    r[2 * len2 - 1] = p;


    for (j = 1; j < 2 * len2; j++)
    {
        a1 = r[j];
        r[j] = r[j] + r[j - 1];
        c = r[j] - a1;
        r[j - 1] = (a1 - (r[j] - c)) + (r[j - 1] - c);
    }

    for (i = 0; i < 2 * len2 - 1; i++)
    {
        res = res + r[i];
    }
    res = res + r[2 * len2 - 1];

    if (err != NULL)
    {
        g = (4 * len2 - 2) * D_EPS / (1 - (4 * len2 - 2) * D_EPS);
        a1 = _d_vec_norm(vec1, len2);
        a2 = _d_vec_norm(vec2, len2);
        *err =
            (D_EPS + 2 * g * g) * fabs(res) + g * g * g * sqrt(a1) * sqrt(a2);
    }

    _d_vec_clear(r);

    return res;
}
