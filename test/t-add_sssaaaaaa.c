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

    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
    int i, j, result;
    flint_rand_t state;
    flint_randinit(state);

    printf("add_sssaaaaaa....");
    fflush(stdout);

    for (i = 0; i < 1000000; i++)
    {
        mp_limb_t s[3], t[3], a[3], b[3];

        for (j = 0; j < 3; j++)
        {
            s[j] = n_randtest(state);
            t[j] = n_randtest(state);
            a[j] = n_randtest(state);
            b[j] = n_randtest(state);
        }

        add_sssaaaaaa(s[2], s[1], s[0], a[2], a[1], a[0], b[2], b[1], b[0]);

        mpn_add_n(t, a, b, 3);

        result = ((s[2] == t[2]) && (s[1] == t[1]) && (s[0] == t[0]));
        if (!result)
        {
            printf("FAIL:\n");
            printf("a[2] = %lu, a[1] = %lu, a[0] = %lu\n", a[2], a[1], a[0]); 
            printf("b[2] = %lu, b[1] = %lu, b[0] = %lu\n", b[2], b[1], b[0]); 
            printf("s[2] = %lu, s[1] = %lu, s[0] = %lu\n", s[2], s[1], s[0]); 
            printf("t[2] = %lu, t[1] = %lu, t[0] = %lu\n", t[2], t[1], t[0]); 
            abort();
        }
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
