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

    Copyright (C) 2015 Vladimir Glazachev
   
******************************************************************************/

#include "aprcl.h"

/*
    Returns the index of divisor p on R factors list.
*/
int
_p_ind(const aprcl_config conf, ulong p)
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
p_power_in_q(ulong q, ulong p)
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

