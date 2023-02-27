/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "aprcl.h"

/*
    Returns the index of divisor p on R factors list.
*/
int
_aprcl_p_ind(const aprcl_config conf, ulong p)
{
    int i;
    for (i = 0; i < conf->rs.num; i++)
        if (p == conf->rs.p[i])
            return i;
    return -1;
}

/*
    Returns k such that p^k | q and p^{k + 1} not | q
*/
ulong
aprcl_p_power_in_q(ulong q, ulong p)
{
    ulong k, q_temp;
    k = 0;
    q_temp = q;
    while (q_temp % p == 0 && q_temp != 0)
    {
        k++;
        q_temp /= p;
    }
    return k;
}
