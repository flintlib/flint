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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "aprcl.h"

int main(void)
{
    ulong i, j;
    FLINT_TEST_INIT(state);
   
    flint_printf("unity_zp_aut_inv....");
    fflush(stdout);

    flint_printf("\n");

    for (i = 0; i < 1; i++)
    {
        ulong p, k, x;
        fmpz_t n;
        unity_zp f, g;

        p = 3;
        k = 1;
        x = 7;

        fmpz_init_set_ui(n, 23);

        unity_zp_init(f, p, k, n);
        unity_zp_init(g, p, k, n);

        jacobi_pq_not2(g, 13, 3);

        unity_zp_aut_inv(f, g, x);
        unity_zp_print(g);
        flint_printf("\n");
        unity_zp_print(f);

        fmpz_clear(n);
        unity_zp_clear(f);
        unity_zp_clear(g);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("NO TEST\n");
    return 0;
}

