/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "ulong_extras.h"
#include "arith.h"

double
arith_dedekind_sum_coprime_d(double h, double k)
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
