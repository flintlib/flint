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
#include "fmpz_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    flint_rand_t state;
    flint_randinit(state);

    printf("set_fmpz_frac....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, y;
        fmpz_t p, q;

        fmpq_init(x);
        fmpq_init(y);
        fmpz_init(p);
        fmpz_init(q);

        fmpz_randtest(p, state, 100);
        fmpz_randtest_not_zero(q, state, 100);

        fmpq_set_fmpz_frac(x, p, q);

        fmpz_set(fmpq_numref(y), p);
        fmpz_set(fmpq_denref(y), q);
        fmpq_canonicalise(y);

        if (!fmpq_is_canonical(x) || !fmpq_equal(x, y))
        {
            printf("FAIL");
            printf("p: "); fmpz_print(p); printf("\n"); 
            printf("q: "); fmpz_print(q); printf("\n"); 
            printf("x: "); fmpq_print(x); printf("\n"); 
            printf("y: "); fmpq_print(y); printf("\n"); 
            abort();
        }

        fmpq_clear(x);
        fmpq_clear(y);
        fmpz_clear(p);
        fmpz_clear(q);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
