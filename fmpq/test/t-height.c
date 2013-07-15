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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"

int
main(void)
{
    int i;
    flint_rand_t state;
    flint_randinit(state);

    printf("height....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        fmpq_t x;
        fmpz_t h;
        mp_bitcnt_t b;

        fmpz_init(h);
        fmpq_init(x);
        fmpq_randtest(x, state, 200);

        fmpq_height(h, x);
        b = fmpq_height_bits(x);

        if (!fmpz_bits(h) == b)
        {
            printf("FAIL!\n");
            printf("x: ");
            fmpq_print(x);
            printf("\nh: ");
            fmpz_print(h);
            printf("\nb: %ld\n", b);
        }

        fmpq_clear(x);
        fmpz_clear(h);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
