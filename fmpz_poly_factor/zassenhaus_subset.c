/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_factor.h"

/* first subset of size m */
void zassenhaus_subset_first(slong * s, slong r, slong m)
{
    slong i;

    FLINT_ASSERT(0 <= m && m <= r);

    for (i = 0; i < r; i++)
    {
        if (i >= m)
            s[i] = (s[i] < 0) ? s[i] : -s[i] - 1;
        else
            s[i] = (s[i] >= 0) ? s[i] : -s[i] - 1;
    }
}

/* next subset of the same size, 1 for success, 0 for failure */
int zassenhaus_subset_next(slong * s, slong r)
{
    slong i, j, k;

    i = 0;
    while (i < r && s[i] < 0)
        i++;
    j = i;

    while (i < r && s[i] >= 0)
        i++;
    k = i;

    if (k == 0 || k >= r)
        return 0;

    s[k] = -s[k] - 1;
    s[k - 1] = -s[k - 1] - 1;

    if (j > 0)
    {
        for (i = 0; i < k - j - 1; i++)
            if (s[i] < 0)
                s[i] = -s[i] - 1;

        for (i = k - j - 1; i < k - 1; i++)
            if (s[i] >= 0)
                s[i] = -s[i] - 1;
    }
    return 1;
}

/* next subset of same size and disjoint from current, delete current idxs */
slong zassenhaus_subset_next_disjoint(slong * s, slong r)
{
    slong i, j, last, total, min;

    total = 0;
    last = r - 1;
    for (i = 0; i < r; i++)
    {
        if (s[i] >= 0)
        {
            total++;
            last = i;
        }
    }

    j = 0;
    for (i = 0; i < r; i++)
        if (s[i] < 0)
            s[j++] = s[i];

    if (r - total < total || total < 1 || last == r - 1)
        return 0;

    min = FLINT_MIN(total - 1, last - total + 1);

    for (i = 0; i < min; i++)
        s[i] = -s[i] - 1;

    for (i = last - total + 1; i < last - min + 1; i++)
        s[i] = -s[i] - 1;

    return 1;
}

