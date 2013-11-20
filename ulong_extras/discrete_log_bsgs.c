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

    Copyright (C) 2013 Mike Hansen
 
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_discrete_log_bsgs(mp_limb_t b, mp_limb_t a, mp_limb_t n)
{
    int i, j;
    double ninv;
    mp_limb_t m, a_m, c;
    mp_limb_t * table;

    ninv = n_precompute_inverse(n);

    m = ceil(sqrt((double) n));

    table = (mp_limb_t*)flint_malloc(m * sizeof(mp_limb_t));

    table[0] = 1;
    for (j = 1; j < m; j++)
    {
        table[j] = n_mulmod_precomp(table[j-1], a, n, ninv);
    }

    a_m = n_invmod(a, n);
    a_m = n_powmod_precomp(a_m, m, n, ninv);

    c = b;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < m; j++)
        {
            if (c == table[j])
            {
                flint_free(table);
                return i * m + j;
            }
        }
        c = n_mulmod_precomp(c, a_m, n, ninv);
    }

    flint_free(table);
    flint_printf("Exception (n_discrete_log_bsgs).  discrete log not found.\n");
    abort();
}
