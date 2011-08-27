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

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "math.h"
#include "arith.h"
#include "ulong_extras.h"

double
dedekind_cosine_sum_d(ulong k, ulong n)
{
    double sum, s, s1, s2;
    double t, u;
    ulong h;
    sum = 0.0;

    if (k <= 2)
    {
        if (k == 0)
            return 0.0;
        else if (k == 1)
            return 1.0;
        else if (n % 2 == 1)
            return -1.0;
        else
            return 1.0;
    }

    t = 3.141592653589793238462643 / (6 * k);

    u = fmod(12.0 * ((double) n), 12 * k);

    for (h = 0; h < (k + 1) / 2; h++)
    {
        if (n_gcd(k, h) == 1UL)
        {
            s = dedekind_sum_coprime_d(h, k);
            s *= (6 * k);
            s = floor(s + 0.5);

            s1 = fmod(s - u*h, 12 * k);
            s2 = fmod(-s - u*(k-h), 12 * k);

            sum += cos(t*s1);
            sum += cos(t*s2);
        }
    }

    return sum;
}
