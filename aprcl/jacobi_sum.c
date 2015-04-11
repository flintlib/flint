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
#include "ulong_extras.h"

void jacobi_pq(unity_root result, ulong q, ulong p)
{
    int i;
    mp_ptr table;

    unity_init(result, n_pow(p, (q - 1) / p));
    table = f_table(q);

    for (i = 0; i < q - 2; i++)
    {
        unity_root temp;
        unity_init(temp, result->power);
        unity_nth_root(temp, i + table[i] + 1);
        unity_roots_add(result, result, temp);
        unity_clear(temp);
    }
}

