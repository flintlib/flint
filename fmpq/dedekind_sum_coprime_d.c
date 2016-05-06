/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "fmpq.h"

double
fmpq_dedekind_sum_coprime_d(double h, double k)
{
    double a, b, t, s, sign;

    if (k <= 2)
        return 0.0;

    a = k;
    b = h;
    s = 0.0;
    sign = 1.0;

    while (b > 0)
    {
        s += sign * (1.0 + a*a + b*b) / (a * b);
        t = fmod(a, b);
        a = b;
        b = t;
        sign = -sign;
    }

    s *= (1./12);

    if (sign < 0)
        s -= 0.25;

    return s;
}

