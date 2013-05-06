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

    printf("cfrac_bound....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, r;
        fmpz * c;
        long n, bound;

        fmpq_init(x);
        fmpq_init(r);

        /* Test worst case (quotient of Fibonacci numbers) */
        if (n_randint(state, 50) == 1)
        {
            long v = 1 + n_randint(state, 1000);
            fmpz_fib_ui(fmpq_numref(x), v + 1);
            fmpz_fib_ui(fmpq_denref(x), v);
        }
        else
        {
            fmpq_randtest(x, state, 1 + n_randint(state, 1000));
        }

        bound = fmpq_cfrac_bound(x);
        c = _fmpz_vec_init(bound + 10);
        n = fmpq_get_cfrac(c, r, x, bound);

        if (n > bound)
        {
            printf("FAIL: length=%ld > bound=%ld\n", n, bound);
            abort();
        }

        _fmpz_vec_clear(c, bound + 10);
        fmpq_clear(x);
        fmpq_clear(r);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
