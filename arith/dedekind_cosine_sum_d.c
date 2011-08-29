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
#include <string.h>
#include <mpir.h>
#include "flint.h"
#include "math.h"
#include "arith.h"
#include "ulong_extras.h"

static __inline__ void
update_count(int * count, double p, double k)
{
    if (p < 0)
        p = -p;

    p = fmod(p, 12.0 * k);

    if (p >= 6.0 * k)
        p = 12.0 * k - p;

    if (p > 3.0 * k)
        count[(int)(6.0*k - p)]--;
    else
        count[(int)(p)]++;
}

double
dedekind_cosine_sum_d(ulong k, ulong n)
{
    double sum, s;
    double t, q;
    int * count;
    ulong i, h;

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

    count = calloc((3*k + 1), sizeof(int));

    t = 3.141592653589793238462643 / (6.0 * k);
    q = fmod(12.0 * ((double) n), 12.0 * k);

    for (h = 0; h < (k + 1) / 2; h++)
    {
        if (n_gcd(k, h) == 1UL)
        {
            s = dedekind_sum_coprime_d(h, k);
            s = floor(s * (6.0 * k) + 0.5);
            update_count(count, s - q*h, k);
            update_count(count, -s - q*(k-h), k);
        }
    }

    sum = count[0];
    for (i = 1; i < 3 * k; i++)
    {
        if (count[i] != 0)
            sum += count[i] * cos(t * i);
    }

    free(count);

    return sum;
}
