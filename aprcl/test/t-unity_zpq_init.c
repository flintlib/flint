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
    int i;
    FLINT_TEST_INIT(state);
   
    flint_printf("unity_zpq_init....");
    fflush(stdout);

    for (i = 0; i < 100; i++)
    {
        ulong p, q;
        fmpz_t n;
        unity_zpq f;

        fmpz_init(n);

        p = n_randprime(state, 2 + n_randint(state, 16), 0);
        q = n_randprime(state, 2 + n_randint(state, 16), 0);

        fmpz_randtest_not_zero(n, state, 200);

        unity_zpq_init(f, q, p, n);
        unity_zpq_clear(f);

        fmpz_clear(n);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

